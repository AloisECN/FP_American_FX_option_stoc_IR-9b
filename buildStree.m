% function [S_tree, jd_S, ju_S, p_hat_S] = buildStree(N, h, kappa, sigma_S, sigma_r, S0, Y0, ro, I, zeta_tree, yearFrac)
% 
%     Y_tree = NaN(N+1, N+1);
%     for i = 1 : N + 1
%         for j = 1 : i
%             Y_tree(i, j) = Y0 + (2*(j - 1) - (i - 1))*sqrt(h);  
%         end
%     end
% 
%     jd_S = NaN(N+1, N+1, N+1);
%     ju_S = NaN(N+1, N+1, N+1);
%     p_hat_S = NaN(N+1, N+1, N+1);
% 
%     for i = 1 : N
%         for j = 1 : i
%             for k = 1 : i
% 
%                 drift = (1 / sqrt(1-ro^2)) * (sigma_r/sigma_S - ro*kappa) * zeta_tree(i, k);
%                 jd_S(i, j, k) = j + floor((1 + sqrt(h) * drift)/2);
%                 jd_S(i, j, k) = max(1, min(i, jd_S(i, j, k)));
%                 ju_S(i, j, k) = jd_S(i, j, k) + 1;
% 
%                 % CORREZIONE: X_tree(i+1, ...) per il nodo futuro
%                 num = drift*h + Y_tree(i, j) - Y_tree(i+1, jd_S(i, j, k));
%                 den = 2*sqrt(h);
%                 p_hat_S(i, j, k) = max(0, min(1, num/den));
%             end
%         end
%     end
% 
%     S_tree = NaN(N+1, N+1, N+1);
%     t_nodes = (0:N) * h;
%     I_ti = interp1(yearFrac, I, t_nodes, 'linear', 'extrap');
% 
%     for i = 1 : N + 1 
%         for j = 1 : i
%             for k = 1 : i
%                 S_tree(i, j, k) = S0 * exp(sigma_S * (sqrt(1-ro^2) * Y_tree(i, j) ...
%                     + ro * zeta_tree(i, k)) + I_ti(i) - (sigma_S^2 / 2) * t_nodes(i));
%             end
%         end
%     end
% 
% end



function [S_tree, jd_S, ju_S, p_hat_S] = buildStree(N, h, kappa, sigma_S, sigma_r, S0, Y0, ro, I, zeta_tree, yearFrac)

    % 1. Vettorializzazione del calcolo di Y_tree (solo ciclo i)
    Y_tree = NaN(N+1, N+1);
    for i = 1 : N + 1
        j_vec = 1:i; % Vettore riga
        Y_tree(i, 1:i) = Y0 + (2*(j_vec - 1) - (i - 1)) * sqrt(h);  
    end
    
    % Inizializzazione output 3D
    jd_S = NaN(N+1, N+1, N+1);
    ju_S = NaN(N+1, N+1, N+1);
    p_hat_S = NaN(N+1, N+1, N+1);
    
    % Precalcolo di costanti per velocizzare il ciclo
    const_drift = (1 / sqrt(1-ro^2)) * (sigma_r/sigma_S - ro*kappa);
    den = 2 * sqrt(h);
    
    % 2. Vettorializzazione delle probabilità e degli indici (cicli j, k eliminati)
    for i = 1 : N
        j_vec = (1:i)'; % Vettore COLONNA (i x 1)
        
        % zeta_tree(i, 1:i) è un vettore RIGA (1 x i). 
        % Il drift dipende solo da k, quindi è un vettore riga (1 x i).
        drift_vec = const_drift * zeta_tree(i, 1:i);
        
        % --- MAGIA DELL'ESPANSIONE IMPLICITA ---
        % Sommando j_vec (i x 1) e drift_vec (1 x i), MATLAB crea in automatico 
        % una matrice (i x i) con tutte le combinazioni incrociate!
        jd_mat = j_vec + floor((1 + sqrt(h) * drift_vec) / 2);
        
        % Applichiamo i limiti a tutta la matrice
        jd_mat = max(1, min(i, jd_mat));
        ju_mat = jd_mat + 1;
        
        % --- TRUCCO DI INDICIZZAZIONE PER IL NODO FUTURO ---
        % Invece di fare Y_tree(i+1, jd_S(i,j,k)), estraiamo l'intera riga i+1...
        Y_next_row = Y_tree(i+1, :);
        % ...e passiamo jd_mat come indice! MATLAB sostituirà ogni indice 
        % in jd_mat con il valore corrispondente presente in Y_next_row.
        Y_next_mat = Y_next_row(jd_mat); % Risultato: matrice (i x i)
        
        Y_curr_vec = Y_tree(i, 1:i)'; % Vettore colonna (i x 1)
        
        % Calcolo probabilità (riga + colonna - matrice = matrice)
        num_mat = drift_vec * h + Y_curr_vec - Y_next_mat;
        p_mat = max(0, min(1, num_mat / den));
        
        % Assegnazione alle "fette" 3D
        jd_S(i, 1:i, 1:i) = jd_mat;
        ju_S(i, 1:i, 1:i) = ju_mat;
        p_hat_S(i, 1:i, 1:i) = p_mat;
    end
    
    % 3. Vettorializzazione del calcolo di S_tree
    S_tree = NaN(N+1, N+1, N+1);
    t_nodes = (0:N) * h;
    I_ti = interp1(yearFrac, I, t_nodes, 'linear', 'extrap');
    
    for i = 1 : N + 1 
        % Y dipende da j -> Vettore COLONNA (i x 1)
        Y_comp = sqrt(1-ro^2) * Y_tree(i, 1:i)';
        
        % zeta dipende da k -> Vettore RIGA (1 x i)
        zeta_comp = ro * zeta_tree(i, 1:i);
        
        % L'espansione implicita crea la matrice (i x i)
        power_mat = sigma_S * (Y_comp + zeta_comp) + I_ti(i) - (sigma_S^2 / 2) * t_nodes(i);
        
        % Calcolo finale e assegnazione alla slice i-esima
        S_tree(i, 1:i, 1:i) = S0 * exp(power_mat);
    end
    
end