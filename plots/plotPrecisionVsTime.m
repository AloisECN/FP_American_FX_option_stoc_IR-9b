function fig = plotPrecisionVsTime(N_vec, prices, times, methodName)

    fig = figure('Name', sprintf('%s: Precision vs Time', methodName), ...
                 'Color', 'w', 'Position', [250, 250, 700, 400]);
    
    green = [0.4660, 0.6740, 0.1880];
    
    yyaxis left;  
    plot(N_vec, prices, '-ok', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); 
    ylabel('ATM Option Price (USD)', 'Interpreter', 'latex', 'FontSize', 12);
    ax = gca;
    ax.YAxis(1).Color = 'k'; 
    
    yyaxis right; 
    plot(N_vec, times, '-^', 'LineWidth', 1.5, 'Color', green, ...
         'MarkerFaceColor', green, 'MarkerEdgeColor', green); 
    ylabel('Computational Time (s)', 'Interpreter', 'latex', 'FontSize', 12);
    ax.YAxis(2).Color = green; 
    
    xlabel('Number of Time Steps ($N$)', 'Interpreter', 'latex', 'FontSize', 12); 
    grid on; grid minor;
    ax.TickLabelInterpreter = 'latex';
    ax.Box = 'on';
    
    title(sprintf('\\textbf{%s: Convergence \\& Computational Cost vs $N$}', methodName), ...
          'Interpreter', 'latex', 'FontSize', 13);    
end