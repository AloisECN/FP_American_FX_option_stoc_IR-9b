function plotCostOfCarry(yearFraction, costOfCarry)
% PLOTCOSTOFCARRY Genera un grafico rapido per controllare la 
% calibrazione del cost of carry.
    
    % Apre una finestra di dimensioni ragionevoli con sfondo bianco
    figure('Name', 'Check Calibrazione Cost of Carry', 'Color', 'w', 'Position', [200, 200, 700, 450]);
    
    % Grafico del Cost of Carry c(t)
    plot(yearFraction, costOfCarry, '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    grid on;
    
    % Impostazioni degli assi e interprete LaTeX
    set(gca, 'TickLabelInterpreter', 'latex');
    xlabel('Time (Years)', 'Interpreter', 'latex', 'FontSize', 11);
    ylabel('$c(t)$', 'Interpreter', 'latex', 'FontSize', 11);
    
    % Titolo del grafico
    title('Instantaneous Cost of Carry $c(t)$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
end