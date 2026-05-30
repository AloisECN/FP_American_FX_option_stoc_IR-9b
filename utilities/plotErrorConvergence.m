function fig = plotErrorConvergence(N_vec, prices, methodName)
    
    if strcmpi(methodName, 'Wei')
        colore = [0, 0.4470, 0.7410]; 
    elseif strcmpi(methodName, 'Robust')
        colore = [0.8500, 0.3250, 0.0980]; 
    else
        colore = [0.4660, 0.6740, 0.1880]; 
    end
    
    % Assume the price with highest N as benchmark
    P_ref = prices(end);
    N_plot = N_vec(1:end-1);
    err_abs = abs(prices(1:end-1) - P_ref);
    
    fig = figure('Name', ['Error - ', methodName], 'Color', 'w', 'Position', [150, 150, 750, 450]);
    
    loglog(N_plot, err_abs, '-s', 'LineWidth', 1.5, ...
        'MarkerSize', 7, 'MarkerFaceColor', colore, 'Color', colore); 
    hold on;
    
    err_ideale_1 = err_abs(1) * (N_plot(1) ./ N_plot);
    loglog(N_plot, err_ideale_1, '--', 'LineWidth', 1.5, 'Color', [0.2, 0.2, 0.2]);
    
    err_ideale_2 = err_abs(1) * (N_plot(1) ./ N_plot).^2;
    loglog(N_plot, err_ideale_2, ':', 'LineWidth', 1.5, 'Color', [0.5, 0.5, 0.5]);
    
    grid on; grid minor;
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    
    title(sprintf('\\textbf{%s Method:} Error Convergence', methodName), ...
        'Interpreter', 'latex', 'FontSize', 13);
    xlabel('Number of Time Steps ($N$)', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Absolute Error $|P_N - P_{max}|$', 'Interpreter', 'latex', 'FontSize', 12);
    
    lgd = legend([methodName, ' Error'], 'Convergence $O(1/N)$', 'Convergence $O(1/N^2)$', ...
           'Location', 'southwest', 'Interpreter', 'latex');
    lgd.Box = 'off';
    
    hold off;
end