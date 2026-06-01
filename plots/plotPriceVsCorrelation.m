function fig = plotPriceVsCorrelation(rho_vec, prices)

    fig = figure('Name', 'Price vs Correlation', 'Color', 'w', 'Position', [200, 200, 700, 400]);
    
    blue = [0, 0.4470, 0.7410];
    
    plot(rho_vec, prices, '-o', 'LineWidth', 1, ...
         'MarkerSize', 3, 'MarkerFaceColor', 'w', ...
         'MarkerEdgeColor', blue, 'Color', blue);
    
    grid on; grid minor;
    ax = gca;
    ax.Box = 'on';
    ax.TickLabelInterpreter = 'latex';
    ax.GridAlpha = 0.2;
    ax.MinorGridAlpha = 0.08;
    
    title('\textbf{American Put Option Price vs Correlation}', 'Interpreter', 'latex', 'FontSize', 13);
    xlabel('Correlation ($\rho$)', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Option Price (USD)', 'Interpreter', 'latex', 'FontSize', 12);   
end