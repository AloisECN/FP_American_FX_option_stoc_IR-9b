function fig = plotPriceVsStrike(K_vec, prices, S0)
    
    colore = [0.8500, 0.3250, 0.0980]; % Robust orange
    
    fig = figure('Name', 'Price vs Strike', 'Color', 'w', 'Position', [150, 150, 750, 450]);
    
    plot(K_vec, prices, '-o', 'LineWidth', 1.5, ...
        'MarkerSize', 7, 'MarkerFaceColor', colore, 'Color', colore);
    hold on;
    
    % Vertical line at S0 (ATM reference)
    xline(S0, '--', 'LineWidth', 1.2, 'Color', [0.2, 0.2, 0.2], ...
        'Label', '$S_0$ (ATM)', 'Interpreter', 'latex', ...
        'LabelVerticalAlignment', 'bottom', 'FontSize', 11);
    
    grid on; grid minor;
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    
    title('\textbf{Robust Method:} Price vs Strike $K$', ...
        'Interpreter', 'latex', 'FontSize', 13);
    xlabel('Strike $K$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Option Price', 'Interpreter', 'latex', 'FontSize', 12);
    
    lgd = legend('Robust Price', 'Location', 'northeast', 'Interpreter', 'latex');
    lgd.Box = 'off';
    
    hold off;
    
end