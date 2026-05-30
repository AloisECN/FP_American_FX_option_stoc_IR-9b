function fig = plotErrorConvergence(N_vec, prices, methodName)
    
    % 1. Scelta automatica del colore
    if strcmpi(methodName, 'Wei')
        colore = [0, 0.4470, 0.7410]; % Blu per Wei
    elseif strcmpi(methodName, 'Robust')
        colore = [0.8500, 0.3250, 0.0980]; % Rosso/Arancio per Robusto
    else
        colore = [0.4660, 0.6740, 0.1880]; % Verde se scrivi altro
    end

    % 2. Calcolo dell'errore (usiamo l'ultimo N come "prezzo vero")
    P_ref = prices(end);
    N_plot = N_vec(1:end-1);
    err_abs = abs(prices(1:end-1) - P_ref);
    
    % 3. Creazione del grafico
    fig = figure('Name', ['Errore - ', methodName], 'Color', 'w', 'Position', [150, 150, 750, 450]);
    
    % I dati empirici
    loglog(N_plot, err_abs, '-s', 'LineWidth', 1.5, ...
        'MarkerSize', 7, 'MarkerFaceColor', colore, 'Color', colore); 
    hold on;
    
    % Retta ideale di riferimento O(1/N) (Tratteggiata nera)
    err_ideale_1 = err_abs(1) * (N_plot(1) ./ N_plot);
    loglog(N_plot, err_ideale_1, '--', 'LineWidth', 1.5, 'Color', [0.2, 0.2, 0.2]);
    
    % NUOVA: Retta ideale di riferimento O(1/N^2) (Punteggiata grigia)
    err_ideale_2 = err_abs(1) * (N_plot(1) ./ N_plot).^2;
    loglog(N_plot, err_ideale_2, ':', 'LineWidth', 1.5, 'Color', [0.5, 0.5, 0.5]);
    
    % 4. Estetica (griglia e testi)
    grid on; grid minor;
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    
    title(sprintf('\\textbf{Metodo %s:} Convergenza dell''Errore', methodName), ...
        'Interpreter', 'latex', 'FontSize', 13);
    xlabel('Numero di passi temporali ($N$)', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Errore Assoluto $|P_N - P_{max}|$', 'Interpreter', 'latex', 'FontSize', 12);
    
    % Legenda aggiornata con le due rette teoriche
    lgd = legend(['Errore ', methodName], 'Convergenza $O(1/N)$', 'Convergenza $O(1/N^2)$', ...
           'Location', 'southwest', 'Interpreter', 'latex');
    lgd.Box = 'off';
    
    hold off;
end