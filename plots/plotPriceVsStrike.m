function fig = plotPriceVsStrike(K_vec, prices, S0)

    fig = figure('Name', 'Price vs Strike', 'Color', 'w', 'Position', [200, 200, 700, 400]);
    
    blue = [0, 0.4470, 0.7410];
    red  = [0.8500, 0.3250, 0.0980];
    
    hold on;
    
    plot(K_vec, prices, '-o', 'LineWidth', 1, ...
         'MarkerSize', 3, 'MarkerFaceColor', 'w', ...
         'MarkerEdgeColor', blue, 'Color', blue);
    
    xline(S0, '--', 'LineWidth', 1.5, 'Color', red); 
    
    grid on; grid minor;
    ax = gca;
    ax.Box = 'on';
    ax.TickLabelInterpreter = 'latex';
    ax.GridAlpha = 0.2;
    ax.MinorGridAlpha = 0.08;
    
    title('\textbf{American Put Option Price vs Strike}', 'Interpreter', 'latex', 'FontSize', 13);
    xlabel('Strike Price ($K$)', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Option Price (USD)', 'Interpreter', 'latex', 'FontSize', 12);
    
    lgd = legend('American Put Price', 'ATM Level ($S_0$)', ...
                 'Location', 'best', 'Interpreter', 'latex');
    lgd.Box = 'off';
    lgd.FontSize = 11;
    
    hold off;
end