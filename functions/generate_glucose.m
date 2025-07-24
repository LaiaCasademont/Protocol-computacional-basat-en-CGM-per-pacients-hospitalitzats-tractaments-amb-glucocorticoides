function glucose_used = generate_glucose(Glucose)
    % --- paràmetres de l’error relatiu en % (igual que abans) ---
    mu    = 16.10700404;
    sigma = 21.54710424;
    max_pct = 2*sigma;   % 2σ de tall

    % filtre exponencial persistent
    persistent last_glucose
    if isempty(last_glucose)
        last_glucose = Glucose;
    end

    % Genera un error truncat (±2σ) i signe aleatori
    while true
        err_pct = normrnd(mu, sigma);
        if abs(err_pct) <= max_pct
            break
        end
    end
    err = (err_pct/100) * (2*randi([0,1]) - 1);

    % applica l’error
    raw = Glucose * (1 + err);

    % no negatius
    raw = max(raw, 0);

    % suavitza amb un EMA molt lleuger
    alpha = 0.2;  % més baix = més suau
    glucose_used = alpha*raw + (1-alpha)*last_glucose;

    % actualitza
    last_glucose = glucose_used;
end
