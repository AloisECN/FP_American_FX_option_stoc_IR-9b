function df_tree = fwdDiscounts_OU(x_tree, kappa, sigma_r, t_nodes, B0_nodes)
% GETLOCALDISCOUNTTREEHW Calcola i fattori di sconto locali esatti B(t_i, t_{i+1})
% per l'intero albero di Hull-White usando la formula chiusa.
%
% INPUTS:
%   x_tree   - Matrice [N+1 x N+1] con i valori del processo OU puro (senza shift alpha)
%   a        - [scalar] mean reversion speed (kappa)
%   sigma    - [scalar] HW volatility (sigma_r)
%   t_nodes  - Vettore [N+1 x 1] dei tempi dell'albero (es. 0, h, 2h...)
%   B0_nodes - Vettore [N+1 x 1] dei discount factor di mercato P(0, t_i)
%
% OUTPUT:
%   df_tree  - Matrice [N+1 x N+1] dei fattori di sconto forward.
%              df_tree(i, k) è lo sconto da applicare al nodo (i,k) per andare al tempo i+1.

    % Inizializza la matrice dei risultati con le stesse dimensioni dell'albero
    [N_rows, N_cols] = size(x_tree);
    N = length(t_nodes) - 1;
    df_tree = NaN(N_rows, N_cols);

    % Ciclo su ogni intervallo temporale [t_i, t_{i+1}]
    for i = 1:N
        t = t_nodes(i);
        T = t_nodes(i+1);
        
        B0_t = B0_nodes(i);
        B0_T = B0_nodes(i+1);

        % 1. Termine B(t,T) deterministico
        B_HW = (1 - exp(-kappa * (T - t))) / kappa;

        % 2. Termine integrale (varianza) dalla tua funzione
        term1 =  1.5;
        term2 = -2   * exp(-kappa*(T-t))   + 0.5 * exp(-2*kappa*(T-t));
        term3 =  2   * exp(-kappa*T)       - 2   * exp(-kappa*t);
        term4 = -0.5 * exp(-2*kappa*T)     + 0.5 * exp(-2*kappa*t);
        integral_term = (sigma_r^2 / kappa^3) * (term1 + term2 + term3 + term4);

        % 3. Estraiamo tutti i nodi attivi al tempo i
        xt = x_tree(i, :);

        % 4. Calcolo vettoriale di B(t,T) per tutti i nodi k contemporaneamente
        FwdB = (B0_T / B0_t) * exp(-xt .* B_HW - 0.5 * integral_term);

        % Salviamo i fattori di sconto nella matrice
        df_tree(i, :) = FwdB;
    end
end