function fig = plotPriceVsCorrelation(rho_vec, prices)
    
    colore = [0.8500, 0.3250, 0.0980]; % Robust orange, consistent with your palette
    
    fig = figure('Name', 'Price vs Correlation', 'Color', 'w', 'Position', [150, 150, 750, 450]);
    
    plot(rho_vec, prices, '-o', 'LineWidth', 1.5, ...
        'MarkerSize', 7, 'MarkerFaceColor', colore, 'Color', colore);
    
    grid on; grid minor;
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    
    title('\textbf{Robust Method:} Price vs Correlation $\rho$', ...
        'Interpreter', 'latex', 'FontSize', 13);
    xlabel('Correlation $\rho$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Option Price', 'Interpreter', 'latex', 'FontSize', 12);
    
    lgd = legend('Robust Price', 'Location', 'best', 'Interpreter', 'latex');
    lgd.Box = 'off';
    
end