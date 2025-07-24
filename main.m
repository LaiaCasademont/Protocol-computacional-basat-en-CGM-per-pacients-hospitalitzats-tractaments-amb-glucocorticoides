% Type 1 Diabetes Simulator using the Hovorka Model
% 
% Description
% This project implements a version of the Hovorka model that allows for 
% 'in silico' testing of both Artificial Pancreas and Multiple Daily 
% Injections therapies. It provides all tools necessary for the analysis 
% of the results. See the documentation file for further information. 
%clear 
% Example
% 1. Go and fill set_simulation.m
% 2. Run main.m
%
% Universitat de Girona - MICElab
% February 2021

%% Purge residual workspace variables 
clear all; clc; close all; 

%% Set simulation
set_simulation;
poolobj = gcp('nocreate');
delete(poolobj);
% Cohort generation (Virtual Patients, Basal, Variability, CF, CR)
if (generate_cohort)
    disp('Generating new virtual patient cohort');
    
    
    %Virtual Patient Generation
    [cohort] = cohort_generation(cohort);
    
    %Basal Tuning for the Virtual Patient Cohort
    cohort.target_fasting_BG = 160;  % mg/dl
    [cohort] = adjust_basal(cohort, simulator, scenario); % Adjust basal for the generated patients

    % CF and CR Tuning for the Virtual Patient Cohort
    [cohort] = adjust_cf(cohort);
   
    [cohort] = adjust_cr(cohort, scenario, simulator);

    % Clear cohort structure and keep only important data TODO
    [cohort] = transform_data(cohort);
    
    % Save all cohorts
    save(path4myOS(strcat('.\virtual_patients\', cohort_name)),'cohort');
else
    load(cohort_path);
end

%% Hardware
hardware.Insulinpump = struct();
hardware.Insulinpump = select_pump(scenario.insulin_pump_model);
hardware.CGM  = struct();
hardware.CGM  = select_CGM(scenario.cgm_model);

%% Scenario 
if (LOAD_SCENARIO)
    load(SCENARIO_PATH);
    simulator.sim_time = scenario.ndays*60*24; % Simulation time [min]
    scenario.exercise.enabled             = false;  % activate exercise or no
    scenario.exercise.management.enabled  = false;  % activate exercise management or no
else    
   [scenario] = generate_scenario(cohort, scenario, simulator, SCENARIO_PATH, hardware);
end

%% Configuration of some variables (history, metrics, hardware)
[hardware, metrics, history] = configuration(simulator, cohort, hardware);

%% Information about the simulation to be carried out
print_simulation_configuration(simulator, scenario, cohort, simulator.controller_name, LOAD_SCENARIO); 

%% Simulation Loop
glucose_production_rate = 0;
tstam = 0:simulator.ts:simulator.sim_time - simulator.ts;
parallel_computing = (length(cohort)>=minimum_npatients_for_parallel)*maximum_nworkers;
% variable copies to allow for parfor and reduce data transfer burden
controller_parfor = hardware.controller;
time_span         = (1:5:simulator.sim_time);
% parfor (subject=1:length(cohort),parallel_computing)
for subject = 1:length(cohort)
    fprintf('Running %s\n',cohort(subject).name);
    % variables defined to allow parfor and reduce data transfer burden
    exercise_management = scenario.exercise.management;
    simulator_parfor  = simulator;
    Ra                = 0;
    Mean_PGUA_1       = 0;
    exercise_end      = 0;
    exercise_start    = 0;
    Tmaxrescue        = 20;
    % Load Virtual patient
    SensitivityFactor = 1; % TODO CHANGE THIS WHEN GENERATING PATIENTS MAYBE?
    patient           = virtualPatient_Hovorka(SensitivityFactor, cohort(subject).params);
    % Controller initialization for each patient
    controller_parfor(subject)          = init_mycontroller(simulator_parfor, scenario, cohort, controller_parfor(subject), subject, patient);
    simulator_parfor.MDI_insulin_type   = controller_parfor(subject).MDI.insulin_type;    
    
    % Initialization of announcement structures
    Meal_announcement     = struct('gramsCHO', 0, 'Ra', 0);
    Exercise_announcement = struct('intensity', 0, 'duration', 0);
    
    % Virtual Patient Initial Conditions
    xguess  = ones(9, 1);
    options = optimset('Diagnostics', 'off', 'Display', 'off');
    [x, fval, exitflag] = fsolve(@(x) model_hovorka_SS(x, cohort(subject).basal*1000/60, patient), xguess, options);
    if (controller_parfor(subject).MDI.enabled)
        [x, fval, exitflag]        = fsolve(@(x) model_hovorka_SS(x, 0, patient), xguess, options);
        long_acting_insulin        = struct;
        long_acting_insulin.adjust = controller_parfor(subject).specific.slow_insulin_dose;
        long_acting_insulin.hour   = controller_parfor(subject).specific.time_injection;
        long_acting_insulin.type   = controller_parfor(subject).MDI.insulin_type;
        x_slowinsulin              = initial_states_slowinsulin(long_acting_insulin,1,scenario.intra_variability);
        Xkm1 = [x', 0, 0, x_slowinsulin', 0, 0, 0, 0,]; % Adding Carbohydrates States + Slow Insulin Compartments
    else
        Xkm1 = [x', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; % No slow acting insulin        
    end
    
    % insulin sensitivity variability
    if (scenario.intra_variability.enabled)
        [pattern_ts]                    = Basal_InsulinVariability_ts(cohort(subject).pattern_insulin_variability_factor_hour,...
                                                    cohort(subject).pattern_basal_insulin_factor_hour, simulator_parfor.sim_time);
        insulin_sensitivity_variability = pattern_ts.insulin_variability_factor ;
    else
        insulin_sensitivity_variability = ones(1, scenario.ndays*288);
    end

    % CGM noise
    if (scenario.commonsensornoise)
        CGM_noise = scenario.CGM_noise;                % noise shared among patients
    else
        CGM_noise = scenario.CGM_noise(subject,:);     % each patient will have a different noise
    end
    
    % meal scenario
    if (scenario.meal.common_scenario)
        meal_Ra        = scenario.meal.Ra(1).Ra_per_patient;
        meal_time_CHO  = scenario.meal.time_CHO(1).time_CHO_per_patient;
        if scenario.meal.misestimation_enable
            meal_cant_CHO  = scenario.meal.cant_CHO_misestimation(1).cant_CHO_misestimation_per_patient;
        else
            meal_cant_CHO  = scenario.meal.cant_CHO(1).cant_CHO_per_patient;
        end 
    else
        meal_Ra        = scenario.meal.Ra(subject).Ra_per_patient;
        meal_time_CHO  = scenario.meal.time_CHO(subject).time_CHO_per_patient;
        if scenario.meal.misestimation_enable
            meal_cant_CHO  = scenario.meal.cant_CHO_misestimation(subject).cant_CHO_misestimation_per_patient;
        else
            meal_cant_CHO  = scenario.meal.cant_CHO(subject).cant_CHO_per_patient;
        end 
    end
    
    % exercise scenario
    exercise_on        = false; % No exercise at the start by default
    exercise_time      = [];
    exercise_intensity = [];
    exercise_duration  = [];
    if (scenario.exercise.enabled)
        if (scenario.exercise.common_scenario && scenario.exercise.enabled)
            exercise_time      = scenario.exercise.time_vector_per_patient(1).time_vector;
            exercise_intensity = scenario.exercise.intensity_vector_per_patient(1).intensity_vector;
            exercise_duration  = scenario.exercise.duration_vector_per_patient(1).duration_vector;
        else
            exercise_time      = scenario.exercise.time_vector_per_patient(subject).time_vector;
            exercise_intensity = scenario.exercise.intensity_vector_per_patient(subject).intensity_vector;
            exercise_duration  = scenario.exercise.duration_vector_per_patient(subject).duration_vector;
        end
    end

    glucose_production_rate=0; 
    for time = tstam 
        sample                = time/simulator_parfor.ts;
        simulator_parfor.time = time;
        Gp                    = Xkm1(7)*MMOL_2_MGDL/patient(17); % Q1/Vg    
        Gp                    = Gp + glucose_production_rate;
        Glucose               = max(min((Gp + CGM_noise(sample+1)),hardware.CGM.sensor_max),hardware.CGM.sensor_min); % noisy CGM measurement 
        if Glucose > 400
            Glucose=400; 
        end

        % Exercise
        if(scenario.exercise.enabled)
            if exercise_management.enabled
                [exercise_time,exercise_management] = Exercise_management(Glucose,time,exercise_time,...
                    meal_time_CHO,exercise_management,simulator,exercise_duration);
            end            
            current_exercise = find(time == exercise_time);
            if ~isempty(current_exercise)
                exercise_on    = true;
                Mean_PGUA_1    = (0.006*(exercise_intensity(current_exercise))^2 + 1.2264*(exercise_intensity(current_exercise)) - 10.1958)*3;
                exercise_start = time;
                exercise_end   = time + exercise_duration(current_exercise);
            end
            if (exercise_on)
                [simulator_parfor, exercise_on] = GetExerciseResponse(time, exercise_start, exercise_end, exercise_on, Mean_PGUA_1, simulator_parfor);
            end
            if (controller_parfor(subject).announce_exercise)
                exercise_in_advance = find((time + controller_parfor(subject).time_in_advance_ann_ex) == exercise_time);
                if ~isempty(exercise_in_advance)
                    Exercise_announcement = struct('intensity',exercise_intensity(exercise_in_advance), 'duration', exercise_duration(exercise_in_advance));
                else
                    Exercise_announcement = struct('intensity', 0, 'duration', 0);
                end
            end
        end
        
        % Meals
        if (scenario.meal.use_mixedmeals)
            Ra                   = meal_Ra(sample + 1);
            input_meal           = 0;
        else
            meal_now = find(time == meal_time_CHO);
            if ~isempty(meal_now)
                Ra         = 0;
                input_meal = meal_cant_CHO(meal_now);
            end
        end
        % glucocorticoids announcemente/adminsitration
        if(scenario.glucocorticoids.enable)
            Glucocorticoids_administration = find(time == scenario.glucocorticoids.intake_time_vector); %This means that in this time instant the glucocorticoids are delivered, and will affect blood glucose and now the controller knows about that!
        else
            Glucocorticoids_administration=0;
            disp('No s.administra');
        end

       if scenario.glucocorticoids.enable

            hour        = mod(floor(time/60), 24);
            intake_time = scenario.glucocorticoids.intake_time;  % 9:00h

            % Definim finestra i pic
            start_effect = intake_time + 4;   % 13:00
            peak_effect  = intake_time + 8;   % 17:00
            end_effect   = intake_time + 14;  % 23:00

            if hour >= start_effect && hour <= end_effect
                dose = scenario.glucocorticoids.dose;  % en mg

                % Efecte mínim més agressiu
                min_x = max(1 - 0.03 * dose, 0.05);

                % Factor d’efecte sobre insulin
                % (manté paràbola inversa per a x_factor)
                relative_time      = (hour - start_effect) / (end_effect - start_effect);
                distance_from_peak = 4 * (relative_time - 0.5).^2;
                x_factor           = min_x + (1 - min_x) * distance_from_peak;

                % Ara la producció de glucosa amb forma triangular:
                % creix linealment de start→peak i disminueix de peak→end
                if hour <= peak_effect
                    frac = (hour - start_effect) / (peak_effect - start_effect);
                else
                    frac = (end_effect - hour) / (end_effect - peak_effect);
                end
                frac = max(min(frac,1),0);  % dins [0,1]

                % Amplitud del peak augmentada
                peak_extra_glucose      = 6 * dose;  
                glucose_production_rate = peak_extra_glucose * frac;

            else
                x_factor = 1.0;
                glucose_production_rate = 0;
            end

        else
            x_factor = 1.0;
            glucose_production_rate = 0;
        end

        % If controller delivers boluses before meals happen
        if (controller_parfor(subject).announce_meals)
            Meal_announcement.Ra = Ra;
            meal_in_advance = find((time + controller_parfor(subject).time_in_advance_meals) == meal_time_CHO);
            if ~isempty(meal_in_advance)
                Meal_announcement.gramsCHO = meal_cant_CHO(meal_in_advance);
            else
                Meal_announcement.gramsCHO = 0;
            end
        end    
        glucose_used=generate_glucose(Glucose);
        % Controller
         [controller_parfor(subject), insulin, meal_bolus, rescueCHO, corrector_bolus, Tmaxrescue, insulin_glucocorticoids, glucose_used] = ...
            run_mycontroller(time, Glucose, glucose_used, Meal_announcement, Exercise_announcement, controller_parfor(subject), Tmaxrescue,Glucocorticoids_administration);
         total_bolus = meal_bolus + corrector_bolus + insulin_glucocorticoids;
        
        % Modify controller values with exercise management parameters
        if scenario.exercise.enabled && exercise_management.enabled
            total_bolus =  total_bolus * (1-exercise_management.bolus_reduction_percentage/100);
            rescueCHO   =  rescueCHO + exercise_management.aditional_carb;                      
            if exercise_management.CSII_basal_reduction && ~simulator.MDI_enabled
                insulin.pump = insulin.pump * (1-exercise_management.CSII_basal_reduction_percentage/100);
                if time == exercise_end
                    exercise_management.CSII_basal_reduction = false;
                end
            end
        end
        
        % Generate inputs to the model 
        [insulin.pump, total_bolus] = pump(total_bolus,insulin.pump,hardware.Insulinpump);           %[insulin.pump(U/h) bolus(U)]
        model_inputs                = [insulin.pump*1000/60 insulin.MDI total_bolus*1000 rescueCHO]; %[mU/min  mU  mU  g]
        model_disturbances          = [Ra simulator_parfor.exercise.M_E_PGU simulator_parfor.exercise.M_E_PIU simulator_parfor.exercise.M_E_HGP];
        MDI_enabled                 = controller_parfor(subject).MDI.enabled;
        
        % Integrate Hovorka Model during Ts minutes
        [~, Xhov] = ode45(@(t,x) model_hovorka(t, x, model_inputs, ...
                    model_disturbances, scenario.intra_variability, ...
                    insulin_sensitivity_variability(sample + 1), patient, ...
                    simulator_parfor,MDI_enabled,Tmaxrescue,x_factor), [time time + simulator_parfor.ts], Xkm1);
        Xkm1 = Xhov(end,:)';
        
        % Save history variables
        history(subject).state(sample + 1, :)            = Xkm1';
        history(subject).Gp(sample + 1)                  = Gp;                                                          % Q1/Vg
        history(subject).Ravmgmin(sample + 1)            = Ra;
        history(subject).TIME(sample + 1)                = time;
        history(subject).CGM(sample + 1)                 = Glucose;                                                %mg/dl
        history(subject).CGBM(sample + 1)                = glucose_used; 
        history(subject).insulin_CSII(sample + 1)        = (insulin.pump/60)*controller_parfor(subject).specific.ts;    % U
        history(subject).meal_bolus(sample + 1)          = meal_bolus;                                                  % U
        history(subject).corrector_bolus(sample + 1)     = corrector_bolus;                                             % U
        history(subject).total_bolus(sample + 1)         = total_bolus;                                                 % U
        history(subject).long_acting_insulin(sample + 1) = insulin.MDI/1000;                                            % U
        history(subject).rescueCHO(sample + 1)           = rescueCHO;                                                   % g
        history(subject).model_inputs(sample + 1,:)      = model_inputs;    
        history(subject).model_disturbances(sample + 1,:)= model_disturbances;   
        history(subject).insulin_glucocorticoids(sample + 1,:)=insulin_glucocorticoids;
    end
    [postprandial_metrics_patient(subject)] = calculate_posprandial_period(meal_time_CHO,...
        history(subject).CGM,history(subject).rescueCHO,history(subject).corrector_bolus);
end
% copy Controller_parfor into hardware struct
hardware.controller = controller_parfor;

%% Outcomes
% computation of metrics per patient taken outside parfor due to dissimilar structures error of unknown origin
for subject = 1:length(cohort)
    % Individual Patient Metrics
    metrics = compute_metrics(subject, metrics, history, hardware, scenario, simulator.ts);
end

% Population metrics
popmetrics = compute_popmetrics(cohort, metrics, history, scenario, hardware.controller, simulator.print_results);
[postprandial_metrics_cohort]   = calculate_posprandial_period_cohort(postprandial_metrics_patient);
popmetrics.postprandial.patient = postprandial_metrics_patient;
popmetrics.postprandial.cohort  = postprandial_metrics_cohort;

% results name
if simulator.MDI_enabled == true
   results_name = strcat(simulator.controller_name,'_',hardware.controller(1).MDI.insulin_type,...
                                    '_',num2str(hardware.controller(1).specific.time_injection),'H');
else
   results_name = strcat(simulator.controller_name);
end
simulator.results_name = strcat(results_name,'.mat');
simulator.results_path = path4myOS(strcat('.\results\', simulator.results_name));

% figures and excel with metrics
exportar_metrics_excel(cohort, metrics, popmetrics, results_name)
outcome_figures;

% Clean the workspace and save data afterwards
clearvars -except cohort hardware history scenario simulator metrics popmetrics results_name;
save(simulator.results_path);

% Save results files to a folder in \results
name_simulation_results = [datestr(datetime, 'dd-mmm-yyyy') '_T' datestr(datetime, 'HHMMSS')];
results_folder_path     = path4myOS(strcat('.\results\', strcat(results_name + "_", name_simulation_results)));
mkdir(results_folder_path);
movefile(path4myOS(".\results\" + simulator.results_name), results_folder_path);
movefile(path4myOS(".\" + strcat(results_name,'.xls')), results_folder_path);
    
% Move AGP reports from \scripts to the results files folder
for subject = 1:size(cohort, 2)
    documentName = strcat('AGP_',regexprep(string(cohort(subject).name),'#','','ignorecase'));
    documentName = strcat(documentName, '.docx');
    movefile(path4myOS(".\scripts\" + documentName), results_folder_path); 
end

% save figures and close
savefig(strcat(results_folder_path,'\CVGA.fig'))
close(gcf)
savefig(strcat(results_folder_path,'\DailyPopulationCGM.fig'))
close(gcf)
savefig(strcat(results_folder_path,'\Glucose_Insulin_CHO.fig'))
close(gcf)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Additional Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation print
function [] = print_simulation_configuration(simulator, scenario, cohort, controller_name, LOAD_SCENARIO)
fprintf(' \n ');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%% Simulator Settings %%%%%%%%%%%%%%%%%%');
disp(strcat('    1. ', 'Total simulation hours: '         , "   "               , num2str(simulator.sim_time/60)));
disp(strcat('    2. ', 'Sampling time (min): '            , "      "            , num2str(simulator.ts)));

disp('%%%%%%%%%%%%%%%%%%% Cohort Settings %%%%%%%%%%%%%%%%%%%%');
disp(strcat('    1. ', 'Total number of VP: '             , "       "            , num2str(size(cohort,2))));
disp(strcat('    2. ', 'Basal Variability: '              , "        "           , 'true')); % TODO check this
disp(strcat('    3. ', 'Intra-Patient Varibility: '       , " "                  , string(scenario.intra_variability.enabled)));

disp('%%%%%%%%%%%%%%%%%% Scenario Settings %%%%%%%%%%%%%%%%%%%');
disp(strcat('    1. ', 'Scenario name: '                   , "            "      , scenario.name));
if(LOAD_SCENARIO) 
    disp(strcat('    2. ', 'Generating new scenario? '     , "  "                , 'no'));
    disp(strcat('    3. ', 'Scenario loaded from: '        , "     "             , 'path'));
else
        disp(strcat('    2. ', 'Generating new scenario? ' , "  "                , 'yes'));
        disp(strcat('    3. ', 'Scenario saved in: '        , "        "         , 'path'));
end
disp(strcat('    4. ', 'Meal Library: '                    , "             "     , 'yes'));
disp(strcat('    5. ', 'Meal contents: '                   , "            "      , mat2str(scenario.meal.daily_contents)));
disp(strcat('    6. ', 'Meal times: '                      , "               "   , mat2str(scenario.meal.daily_intake_time)));
if (scenario.exercise.enabled)
    disp(strcat('    7. ', 'Exercise: '                    , "                 " , 'yes'));
    disp(strcat('    8. ', 'Exercise days: '               , "            "      , mat2str(scenario.exercise.days)));
    disp(strcat('    9. ', 'Exercise times: '              , "           "       , mat2str(scenario.exercise.days_hour)));
else   
    disp(strcat('    7. ', 'Exercise: '                    , "                 " , 'No'));
end

disp('%%%%%%%%%%%%%%%%%% Hardware Settings %%%%%%%%%%%%%%%%%%%');
disp(strcat('    1. ', 'Controller: '  , " "       , controller_name)); 
disp(strcat('    2. ', 'CGM: '         , "        ", scenario.cgm_model));
disp(strcat('    3. ', 'Pump: '        , "       " , scenario.insulin_pump_model));
if(simulator.MDI_enabled)
    disp(strcat('    4. ', 'Therapy: ' , "    "    , 'MDI'));
else
    disp(strcat('    4. ', 'Therapy: ' , "    "    , 'CSII'));
end
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

fprintf(' \n ');
disp('Press Enter to continue or Ctrl + C to stop the simulation.');
commandwindow;
pause;
end

%% Exercise Model GetExerciseResponse(i, exercise_start, exercise_end, exercise_on, Mean_PGUA_1, scenario, simulator)
function [simulator, exercise_on] = GetExerciseResponse(i, exercise_start, exercise_end, exercise_on, Mean_PGUA_1, simulator)
CURRENT_TIME      = i;
POSTEXERCISE_TIME = 5*60;

if CURRENT_TIME >= exercise_start && CURRENT_TIME <= exercise_end
    syms yy_1(t); YY_1 = dsolve(diff(yy_1) == -(1/30)*yy_1 + (1/30)*Mean_PGUA_1, yy_1(0) == 0.1);
    t = (CURRENT_TIME - exercise_start);
    simulator.exercise.PGUA_1_Act = subs(YY_1); simulator.exercise.PGUA_1_Act = double(simulator.exercise.PGUA_1_Act);
    simulator.exercise.M_E_PGU = 1 + simulator.exercise.PGUA_1_Act*simulator.exercise.PAMM/35;
    simulator.exercise.M_E_PIU = 1 + 2.4*simulator.exercise.PAMM;
    simulator.exercise.M_E_HGP = 1 + simulator.exercise.PGUA_1_Act*simulator.exercise.PAMM/155;
    
elseif CURRENT_TIME > exercise_end && CURRENT_TIME < exercise_end + POSTEXERCISE_TIME 
    syms zz_1(t); ZZ_1 = dsolve(diff(zz_1) == -(1/30)*zz_1 , zz_1(0) == simulator.exercise.PGUA_1_Act);
    t = (CURRENT_TIME - exercise_start);
    simulator.exercise.PGUA_1_NoAct = subs(ZZ_1); simulator.exercise.PGUA_1_NoAct = double(simulator.exercise.PGUA_1_NoAct);
    simulator.exercise.M_E_PGU = 1 + simulator.exercise.PGUA_1_NoAct*simulator.exercise.PAMM/35;
    simulator.exercise.M_E_PIU = 1;
    simulator.exercise.M_E_HGP = 1 + simulator.exercise.PGUA_1_NoAct*simulator.exercise.PAMM/155;
    
else
    simulator.exercise.M_E_PGU  = 1;
    simulator.exercise.M_E_PIU  = 1;
    simulator.exercise.M_E_HGP  = 1;
    exercise_on = false;
end
end

%% Compute metrics
function [metrics] = compute_metrics(subject, metrics, history, hardware, scenario, ts)
% Glycemic Metrics
metrics(subject).mean_CGM         = mean(history(subject).CGM);
metrics(subject).STD_CGM          = std(history(subject).CGM);
metrics(subject).median_CGM       = median(history(subject).CGM);
metrics(subject).CV_CGM           = (metrics(subject).STD_CGM / metrics(subject).mean_CGM)*100;
metrics(subject).prctile25_75_CGM = prctile([history(subject).CGM], [25 75]);
metrics(subject).min_CGM          = min(history(subject).CGM);
metrics(subject).max_CGM          = max(history(subject).CGM);

% Insulin
metrics(subject).Basal_totalinsulindelivered       = sum(history(subject).insulin_CSII) + sum(history(subject).long_acting_insulin); % [U]
metrics(subject).Basal_totalinsulindeliveredperday = metrics(subject).Basal_totalinsulindelivered/scenario.ndays; % [U/day]

metrics(subject).Bolus_totalinsulindelivered       = sum(history(subject).total_bolus); % [U]
metrics(subject).Bolus_totalinsulindeliveredperday = metrics(subject).Bolus_totalinsulindelivered/scenario.ndays; % [U/day]

metrics(subject).totalinsulindelivered             = metrics(subject).Basal_totalinsulindelivered + metrics(subject).Bolus_totalinsulindelivered; % [U]
metrics(subject).totalinsulindeliveredperday       = metrics(subject).totalinsulindelivered/scenario.ndays ; % [U/day]
metrics(subject).STD_totalinsulindeliveredperday   = std(sum((reshape((history(subject).insulin_CSII + history(subject).total_bolus),288,[]))));

% Carbohydrates
metrics(subject).totalgramsrescuecarbs        = sum(history(subject).rescueCHO);
metrics(subject).totaleventsrescuecarbs       = length(find(history(subject).rescueCHO ~= 0));
metrics(subject).totalgramsrescuecarbsperday  = metrics(subject).totalgramsrescuecarbs/scenario.ndays;
metrics(subject).totaleventsrescuecarbsperday = metrics(subject).totaleventsrescuecarbs/scenario.ndays;

% Risk Index
[metrics(subject).RI_CGM, metrics(subject).LBGI_CGM, metrics(subject).HBGI_CGM] = ...
    risk_index(history(subject).CGM);

% Time in Range
[metrics(subject).CGMIN70_140, metrics(subject).CGMIN70_180, metrics(subject).CGMABOVE250, ...
 metrics(subject).CGMIN180_250, metrics(subject).CGMIN70_54, metrics(subject).CGMUNDER54] = ...
 time_in_target(history(subject).CGM);  

[metrics(subject).STD_CGMIN70_180, metrics(subject).STD_CGMABOVE250, ...
 metrics(subject).STD_CGMIN180_250, metrics(subject).STD_CGMIN70_54, metrics(subject).STD_CGMUNDER54] = ...
 STD_time_in_range(history(subject).CGM,ts);  

% Hypoglycemia Events
[metrics(subject).CGMHypos, metrics(subject).CGMTimeHypos] = ...
    hypos_counter(history(subject).CGM, 70, 15, false);

metrics(subject).mean_mealperday  = mean(scenario.meal.meal_per_day(subject).meal_per_day_per_patient);
metrics(subject).STD_mealperday   = std(scenario.meal.meal_per_day(subject).meal_per_day_per_patient);
end

function [hardware] = select_CGM(CGM_model)
    % adapted from UVa simulator v32
    
    if (strcmp(CGM_model,'default'))

        hardware.sensor_PACF=0.7;
        hardware.sensor_type='SU';
        hardware.sensor_gamma=-0.5444;
        hardware.sensor_lambda=15.9574/1.5;   
        hardware.sensor_delta=1.6898;
        hardware.sensor_xi=-5.47;
        hardware.sensor_sampling=1;
        hardware.sensor_min=32;
        hardware.sensor_max=600;
        
    elseif (strcmp(CGM_model,'dexcom'))

        hardware.sensor_PACF=0.7;
        hardware.sensor_type='SU';
        hardware.sensor_gamma=-0.5444;
        hardware.sensor_lambda=15.9574;
        hardware.sensor_delta=1.6898;
        hardware.sensor_xi=-5.47;
        hardware.sensor_sampling=3;
        hardware.sensor_min=39;
        hardware.sensor_max=600; 

    elseif (strcmp(CGM_model,'dexcom25'))

        hardware.sensor_PACF=0.25;
        hardware.sensor_type='SU';
        hardware.sensor_gamma=-0.5444;
        hardware.sensor_lambda=15.9574;
        hardware.sensor_delta=1.6898;
        hardware.sensor_xi=-5.47;
        hardware.sensor_sampling=3;
        hardware.sensor_min=39;
        hardware.sensor_max=600; 
        
    elseif (strcmp(CGM_model,'dexcom50'))

        hardware.sensor_PACF=0.50;
        hardware.sensor_type='SU';
        hardware.sensor_gamma=-0.5444;
        hardware.sensor_lambda=15.9574;
        hardware.sensor_delta=1.6898;
        hardware.sensor_xi=-5.47;
        hardware.sensor_sampling=3;
        hardware.sensor_min=39;
        hardware.sensor_max=600; 
        
     elseif (strcmp(CGM_model,'dexcom70'))

        hardware.sensor_PACF=0.50;
        hardware.sensor_type='SU';
        hardware.sensor_gamma=-0.5444;
        hardware.sensor_lambda=15.9574;
        hardware.sensor_delta=1.6898;
        hardware.sensor_xi=-5.47;
        hardware.sensor_sampling=3;
        hardware.sensor_min=39;
        hardware.sensor_max=600; 
        
     elseif (strcmp(CGM_model,'guardianRT'))
        
        hardware.sensor_PACF=0.7;
        hardware.sensor_type='SU';
        hardware.sensor_gamma=-0.5444;
        hardware.sensor_lambda=15.9574;
        hardware.sensor_delta=1.6898;
        hardware.sensor_xi=-5.47;
        hardware.sensor_sampling=5;
        hardware.sensor_min=39;
        hardware.sensor_max=600;
        
     elseif (strcmp(CGM_model,'navigator'))
        
        hardware.sensor_PACF=0.7;
        hardware.sensor_type='SU';
        hardware.sensor_gamma=-0.5444;
        hardware.sensor_lambda=15.9574;
        hardware.sensor_delta=1.6898;
        hardware.sensor_xi=-5.47;
        hardware.sensor_sampling=1;
        hardware.sensor_min=32;
        hardware.sensor_max=600;

    elseif (strcmp(CGM_model,'navigator25'))
        
        hardware.sensor_PACF=0.25;
        hardware.sensor_type='SU';
        hardware.sensor_gamma=-0.5444;
        hardware.sensor_lambda=15.9574;
        hardware.sensor_delta=1.6898;
        hardware.sensor_xi=-5.47;
        hardware.sensor_sampling=1;
        hardware.sensor_min=32;
        hardware.sensor_max=600;
        
    elseif (strcmp(CGM_model,'navigator50'))
        
        hardware.sensor_PACF=0.5;
        hardware.sensor_type='SU';
        hardware.sensor_gamma=-0.5444;
        hardware.sensor_lambda=15.9574;
        hardware.sensor_delta=1.6898;
        hardware.sensor_xi=-5.47;
        hardware.sensor_sampling=1;
        hardware.sensor_min=32;
        hardware.sensor_max=600;
        
    elseif (strcmp(CGM_model,'navigator70'))
        
        hardware.sensor_PACF=0.7;
        hardware.sensor_type='SU';
        hardware.sensor_gamma=-0.5444;
        hardware.sensor_lambda=15.9574;
        hardware.sensor_delta=1.6898;
        hardware.sensor_xi=-5.47;
        hardware.sensor_sampling=1;
        hardware.sensor_min=32;
        hardware.sensor_max=600;
    
    else
        disp('ERROR: CGM model not defined');
    end

end

function hardware = select_pump(pump_model)
    % adapted from UVa simulator v32
    
    if (strcmp(pump_model,'default'))

        hardware.pump_noise=0;
        hardware.pump_char=1;

        hardware.minbolus=0;        %U
        hardware.maxbolus=75;       %U
        hardware.incbolus=0.05;     %U
        hardware.minbasal=0;        %U/h
        hardware.maxbasal=35;       %U/h
        hardware.incbasal=0.05;     %U/h
        hardware.sampling=1;
        hardware.accuracy_bolus_amount=[0 100];
        hardware.bolus_mean = [0 100];
        hardware.bolus_std2r= [0 0];
        
    elseif (strcmp(pump_model,'cozmo'))
        
        hardware.pump_noise=0;
        hardware.pump_char=1;

        hardware.minbolus=0;
        hardware.maxbolus=75;
        hardware.incbolus=0.05;
        hardware.minbasal=0;
        hardware.maxbasal=35;
        hardware.incbasal=0.05;
        hardware.sampling=1;
        hardware.accuracy_bolus_amount=[0 100];
        hardware.bolus_mean = [0 100];
        hardware.bolus_std2r= [0 0];
        
    elseif (strcmp(pump_model,'insulet'))
        
        hardware.pump_noise=1;
        hardware.pump_char=1;

        hardware.minbolus=0;
        hardware.maxbolus=30;
        hardware.incbolus=0.05;
        hardware.minbasal=0;
        hardware.maxbasal=30;
        hardware.incbasal=0.05;
        hardware.sampling=1;
        hardware.accuracy_bolus_amount=[0 0.05 0.1 0.2 1 35];
        hardware.bolus_mean = [0 0.04965 0.9991 0.2 1.001 35];
        hardware.bolus_std2r= [0.12 0.12 0.08 0.0633  0.04  0];
        

    elseif (strcmp(pump_model,'generic'))
        
        hardware.pump_noise=0;
        hardware.pump_char=1;

        hardware.minbolus=0;
        hardware.maxbolus=75;
        hardware.incbolus=0.05;
        hardware.minbasal=0;
        hardware.maxbasal=35;
        hardware.incbasal=0.05;
        hardware.sampling=1;
        hardware.accuracy_bolus_amount=[0 100];
        hardware.bolus_mean = [0 100];
        hardware.bolus_std2r= [0 0];

    elseif (strcmp(pump_model,'minimal_error'))
        
        hardware.pump_noise=0;
        hardware.pump_char=1;

        hardware.minbolus=0;
        hardware.maxbolus=Inf;
        hardware.incbolus=0.0001;
        hardware.minbasal=0;
        hardware.maxbasal=Inf;
        hardware.incbasal=0.0001;
        hardware.sampling=1;
        hardware.accuracy_bolus_amount=[0 100];
        hardware.bolus_mean = [0 100];
        hardware.bolus_std2r= [0 0];


    elseif (strcmp(pump_model,'no_error'))
        
        hardware.pump_noise=0;
        hardware.pump_char=0;

        hardware.minbolus=0;
        hardware.maxbolus=Inf;
        hardware.incbolus=0.0001;
        hardware.minbasal=0;
        hardware.maxbasal=Inf;
        hardware.incbasal=0.0001;
        hardware.sampling=1;
        hardware.accuracy_bolus_amount=[0 100];
        hardware.bolus_mean = [0 100];
        hardware.bolus_std2r= [0 0];
        
    else
        disp('ERROR: pump model not defined');
    end

end

function [basal_pumped, bolus_pumped]  = pump(bolus,basal,hardware)
    % pump saturation, quantization and random error
    % adapted from UVa simulator v32
    % basal in U/min and bolus in U
       
    if (hardware.pump_char)
      
        % saturate
        bolus_sat = min(max(hardware.minbolus,bolus),hardware.maxbolus);
        basal_sat = min(max(hardware.minbasal,basal),hardware.maxbasal);
        
        % quantize
        qbolus       = hardware.incbolus;
        bolus_pumped = qbolus*floor(bolus_sat/qbolus);
        qbasal       = hardware.incbasal;
        basal_pumped = qbasal*floor(basal_sat/qbasal);
    else
        bolus_pumped = bolus;
        basal_pumped = basal;
    end
   
end

function [STD_IN70_180, STD_ABOVE250, STD_IN180_250, STD_IN70_54, STD_UNDER54]=STD_time_in_range(G,ts)
period              = 24;           %H
samples_per_period  = period*60/ts;     
cant_filas          = length(G)/samples_per_period;
glucose_per_day     = transpose(reshape(G,samples_per_period,[]));
    for dia=1:1:cant_filas
        ABOVE250 = (numel(find(glucose_per_day (dia,:)>250))/ samples_per_period)*100;
        ABOVE180 = (numel(find(glucose_per_day(dia,:)>180))/samples_per_period)*100;
        ABOVE140 = (numel(find(glucose_per_day(dia,:)>140))/samples_per_period)*100;
        UNDER60  = (numel(find(glucose_per_day(dia,:)<60))/samples_per_period)*100;
        UNDER70  = (numel(find(glucose_per_day(dia,:)<70))/samples_per_period)*100;
        UNDER54  = (numel(find(glucose_per_day(dia,:)<54))/samples_per_period)*100;
        
        IN70_180 = 100-(ABOVE180 + UNDER70);
        IN180_250 = ABOVE180 - ABOVE250;
        IN70_140  = 100 - (ABOVE140 + UNDER70);
        IN60_70   = UNDER70 - UNDER60;
        IN54_70   = UNDER70 - UNDER54;
        
        %rangos=[rangos UNDER54, IN54_70, IN70_180, IN180_250, ABOVE250];
        time_ranges_percent_per_day(dia,:) = [UNDER54, IN54_70, IN70_180, IN180_250, ABOVE250];
    end
STD_UNDER54   = std(time_ranges_percent_per_day(:,1));
STD_IN70_54   = std(time_ranges_percent_per_day(:,2));
STD_IN70_180  = std(time_ranges_percent_per_day(:,3));
STD_IN180_250 = std(time_ranges_percent_per_day(:,4));
STD_ABOVE250  = std(time_ranges_percent_per_day(:,5));

end

function [postprandial_metrics_patient] = calculate_posprandial_period(time_meal,CGM,rescue, corrector_bolus)
ts=5;
meal_ts_instances = time_meal./ts+1;
postprandial_metrics_patient.total_CGM_postprandial_vector = [];
postprandial_metrics_patient.total_postprandial_TIR        = [];
for meal_index=1: length(meal_ts_instances)
    if meal_index == length(meal_ts_instances)
        samples_betwen_meals = length(CGM)-meal_ts_instances(end);
    else
        samples_betwen_meals = diff(meal_ts_instances(meal_index:meal_index+1));
    end
    [postprandial_metrics_patient.postprandial_metrics_per_meal(meal_index)] = calculate_posprandial_period_per_meal(meal_ts_instances(meal_index),...
        samples_betwen_meals,CGM,rescue,corrector_bolus);
    postprandial_metrics_patient.total_CGM_postprandial_vector = [postprandial_metrics_patient.total_CGM_postprandial_vector;...
        postprandial_metrics_patient.postprandial_metrics_per_meal(meal_index).CGM_postprandial_vector];
end
[IN70_140, IN70_180, ABOVE250, IN180_250, IN70_54, UNDER54] = time_in_target(postprandial_metrics_patient.total_CGM_postprandial_vector);
postprandial_metrics_patient.total_postprandial_TIR   = [ABOVE250, IN180_250,IN70_180, IN70_54, UNDER54];
end

function  [postprandial_metrics_cohort] = calculate_posprandial_period_cohort(postprandial_metrics_patient)
postprandial_metrics_cohort.total_CGM_postprandial_vector  = [];
for patient=1: length(postprandial_metrics_patient)
    postprandial_metrics_cohort.total_CGM_postprandial_vector  = [postprandial_metrics_cohort.total_CGM_postprandial_vector...
        ; postprandial_metrics_patient(patient).total_CGM_postprandial_vector];
end
[IN70_140, IN70_180, ABOVE250, IN180_250, IN70_54, UNDER54] = time_in_target(postprandial_metrics_cohort.total_CGM_postprandial_vector);
postprandial_metrics_cohort.total_postprandial_TIR   = [ABOVE250, IN180_250,IN70_180, IN70_54, UNDER54];
end

function [meal] = calculate_posprandial_period_per_meal(meal_ts_instances,samples_betwen_meals,CGM,rescue,corrector_bolus)
samples_postprandial_period_4H = 4*60/5;  % 4 h de periodo postprandial
if samples_betwen_meals < samples_postprandial_period_4H
    meal.samples_postprandial_period = samples_betwen_meals;
else
    meal.samples_postprandial_period = samples_postprandial_period_4H;
end
vector_samples = (meal_ts_instances+1):1:meal_ts_instances + meal.samples_postprandial_period;

meal.CGM_postprandial_vector               = CGM(vector_samples);
meal.amount_postprandial_carb_rescue       = length(find(rescue(vector_samples)>0));
meal.minutes_between_meal_carb             = diff([meal_ts_instances meal_ts_instances+find(rescue(vector_samples)>0)'])*5;
meal.amount_postprandial_corrector_bolus   = length(find(corrector_bolus(vector_samples)>0));
meal.minutes_between_meal_corrector_bolus  = diff([meal_ts_instances meal_ts_instances+find(corrector_bolus(vector_samples)>0)'])*5;
[IN70_140, IN70_180, ABOVE250, IN180_250, IN70_54, UNDER54] = time_in_target(meal.CGM_postprandial_vector);
meal.postprandial_TIR   = [ABOVE250, IN180_250,IN70_180, IN70_54, UNDER54];
end

function [hardware, metrics, history] = configuration(simulator, cohort, hardware)
%% Configure history
sim_samples = numel(0:simulator.ts:simulator.sim_time - simulator.ts);
for subject = 1:size(cohort, 2)
    history(subject).name                 = cohort(subject).name;
    history(subject).state                = zeros(sim_samples, 20); % Full state
    history(subject).Gp                   = zeros(sim_samples, 1);  % Plasma glucose
    history(subject).Ravmgmin             = zeros(sim_samples, 1);  % Rate of absorption (mg/min)
    history(subject).TIME                 = zeros(sim_samples, 1);  % Time (min)
    history(subject).CGM                  = zeros(sim_samples, 1);  % [mg/dl]
    history(subject).CGBM                 = zeros(sim_samples, 1); 
    history(subject).insulin_CSII         = zeros(sim_samples, 1);  % U/h
    history(subject).long_acting_insulin  = zeros(sim_samples, 1);  % Glargine infusion if MDI enabled
    history(subject).corrector_bolus      = zeros(sim_samples, 1);  % U/h
    history(subject).meal_bolus           = zeros(sim_samples, 1);  % U/h
    history(subject).total_bolus          = zeros(sim_samples, 1);  % U/h
    history(subject).rescueCHO            = zeros(sim_samples, 1);  % rescue carbs (mg)
    history(subject).insulin_glucocorticoids            = zeros(sim_samples, 1);  % 
end

%% Configure hardware
nPatients = size(cohort, 2);
hardware.controller = struct('name'                  , cell(1, nPatients),...
    'controllerpath'        , cell(1, nPatients),...
    'announce_meals'        , cell(1, nPatients),...
    'time_in_advance_meals' , cell(1, nPatients),...
    'announce_exercise'     , cell(1, nPatients),...
    'time_in_advance_ann_ex', cell(1, nPatients),...
    'MDI'                   , cell(1, nPatients),...
    'history'               , cell(1, nPatients),...
    'specific'              , cell(1, nPatients));

%% Configure metrics
metrics = struct('mean_CGM'                                          , num2cell(zeros(1, nPatients)),...
    'STD_CGM'                                           , num2cell(zeros(1, nPatients)),...
    'median_CGM'                                        , num2cell(zeros(1, nPatients)),...
    'prctile25_75_CGM'                                  , num2cell(zeros(1, nPatients)),...
    'min_CGM'                                           , num2cell(zeros(1, nPatients)),...
    'max_CGM'                                           , num2cell(zeros(1, nPatients)),...
    'CV_CGM'                                            , num2cell(zeros(1, nPatients)),...
    'Basal_totalinsulindelivered'                       , num2cell(zeros(1, nPatients)),...
    'Basal_totalinsulindeliveredperday'                 , num2cell(zeros(1, nPatients)),...
    'Bolus_totalinsulindelivered'                       , num2cell(zeros(1, nPatients)),...
    'Bolus_totalinsulindeliveredperday'                 , num2cell(zeros(1, nPatients)),...
    'totalinsulindelivered'                             , num2cell(zeros(1, nPatients)),...
    'totalinsulindeliveredperday'                       , num2cell(zeros(1, nPatients)),...
    'STD_totalinsulindeliveredperday'                   , num2cell(zeros(1, nPatients)),...
    'totalgramsrescuecarbs'                             , num2cell(zeros(1, nPatients)),...
    'totaleventsrescuecarbs'                            , num2cell(zeros(1, nPatients)),...
    'totalgramsrescuecarbsperday'                       , num2cell(zeros(1, nPatients)),...
    'totaleventsrescuecarbsperday'                      , num2cell(zeros(1, nPatients)),...
    'mean_mealperday'                                   , num2cell(zeros(1, nPatients)),...
    'std_mealperday'                                    , num2cell(zeros(1, nPatients)),...
    'RI_CGM'                                            , num2cell(zeros(1, nPatients)),...
    'LBGI_CGM'                                          , num2cell(zeros(1, nPatients)),...
    'HBGI_CGM'                                          , num2cell(zeros(1, nPatients)),...
    'CGMIN70_140'                                       , num2cell(zeros(1, nPatients)),...
    'CGMIN70_180'                                       , num2cell(zeros(1, nPatients)),...
    'STD_CGMIN70_180'                                   , num2cell(zeros(1, nPatients)),...
    'CGMABOVE250'                                       , num2cell(zeros(1, nPatients)),...
    'STD_CGMABOVE250'                                   , num2cell(zeros(1, nPatients)),...
    'CGMIN180_250'                                      , num2cell(zeros(1, nPatients)),...
    'STD_CGMIN180_250'                                  , num2cell(zeros(1, nPatients)),...
    'CGMIN70_54'                                        , num2cell(zeros(1, nPatients)),...
    'STD_CGMIN70_54'                                    , num2cell(zeros(1, nPatients)),...
    'CGMUNDER54'                                        , num2cell(zeros(1, nPatients)),...
    'STD_CGMUNDER54'                                    , num2cell(zeros(1, nPatients)),...
    'CGMHypos'                                          , num2cell(zeros(1, nPatients)),...
    'CGMTimeHypos'                                      , num2cell(zeros(1, nPatients)));
end