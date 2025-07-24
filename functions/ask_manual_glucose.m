function glucose_used = ask_manual_glucose(Glucose)
    while true
        CGBM_input = input('Introdueix el valor manual de glucosa (mg/dL): ', 's');
        CGBM = str2double(CGBM_input);
        if isempty(CGBM_input) || isnan(CGBM)
            disp('Valor no vàlid. Torna-ho a intentar.');
        else
            % Un cop tenim un valor vàlid, calculem la diferència i decidim
            diff_percent = abs(CGBM - Glucose) / Glucose;
            if diff_percent > 0.25
                disp(['La glucosa CGM= ', num2str(Glucose), '.']);
                disp(['La glucosa CGBM= ', num2str(CGBM_input), '.']);
                disp('Diferència > 25%. S’utilitzarà el valor manual de glucosa.');
                glucose_used = CGBM;
            else
                disp(['La glucosa CGM= ', num2str(Glucose), '.']);
                disp(['La glucosa CGBM= ', num2str(CGBM_input), '.']);
                disp('Diferència <= 25%. S’utilitzarà el valor automàtic de glucosa.');
                glucose_used = Glucose;
            end
            break; % sortim del bucle després d’aquesta única comprovació
        end
    end
end
