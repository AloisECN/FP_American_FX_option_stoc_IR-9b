function [S_tree, jd_S, ju_S, p_hat_S] = buildStree(N, h, kappa, sigma_S, sigma_r, S0, Y0, ro, I, zeta_tree, yearFrac)

    Y_tree = NaN(N+1, N+1);
    for i = 1 : N + 1
        for j = 1 : i
            Y_tree(i, j) = Y0 + (2*(j - 1) - (i - 1))*sqrt(h);  
        end
    end
    
    jd_S = NaN(N+1, N+1, N+1);
    ju_S = NaN(N+1, N+1, N+1);
    p_hat_S = NaN(N+1, N+1, N+1);
    
    for i = 1 : N
        for j = 1 : i
            for k = 1 : i
            
                drift = (1 / sqrt(1-ro^2)) * (sigma_r/sigma_S - ro*kappa) * zeta_tree(i, k);
                jd_S(i, j, k) = j + floor((1 + sqrt(h) * drift)/2);
                jd_S(i, j, k) = max(1, min(i, jd_S(i, j, k)));
                ju_S(i, j, k) = jd_S(i, j, k) + 1;
            
                % CORREZIONE: X_tree(i+1, ...) per il nodo futuro
                num = drift*h + Y_tree(i, j) - Y_tree(i+1, jd_S(i, j, k));
                den = 2*sqrt(h);
                p_hat_S(i, j, k) = max(0, min(1, num/den));
            end
        end
    end
    
    S_tree = NaN(N+1, N+1, N+1);
    t_nodes = (0:N) * h;
    I_ti = interp1(yearFrac, I, t_nodes, 'linear', 'extrap');
    
    for i = 1 : N + 1 
        for j = 1 : i
            for k = 1 : i
                S_tree(i, j, k) = S0 * exp(sigma_S * (sqrt(1-ro^2) * Y_tree(i, j) ...
                    + ro * zeta_tree(i, k)) + I_ti(i) - (sigma_S^2 / 2) * t_nodes(i));
            end
        end
    end
    
end



%% VECTORIALIZED FUNCTION -> MORE EFFICIENT

% function [S_tree, jd_S, ju_S, p_hat_S] = buildStree(N, h, kappa, sigma_S, sigma_r, S0, X0, ro, I, zeta_tree)
% 
%     
%     tempi_idx  = (0:N)'; 
%     spazio_idx = (0:N);  
% 
%     X_tree = X0 + (2 * spazio_idx - tempi_idx) * sqrt(h);
%     X_tree(spazio_idx > tempi_idx) = NaN; % Elimina i nodi impossibili
% 
%    
%     jd_S    = NaN(N+1, N+1, N+1);
%     ju_S    = NaN(N+1, N+1, N+1);
%     p_hat_S = NaN(N+1, N+1, N+1);
% 
%     drift_costante = (1 / sqrt(1 - ro^2)) * (sigma_r / sigma_S - ro * kappa);
% 
%     for i = 0:N-1
% 
%         nodi_S = (0:i)'; 
%         nodi_r = (0:i);  
% 
%         drift_corrente = drift_costante * zeta_tree(i+1, 1:i+1); 
% 
%         salto_giu = nodi_S + floor((1 + sqrt(h) * drift_corrente) / 2);
%         salto_giu = max(0, min(i, salto_giu));
%         salto_su  = salto_giu + 1;
% 
%         X_corrente  = X_tree(i+1, 1:i+1)'; 
%         X_futuro_giu = X0 + (2 * salto_giu - (i + 1)) * sqrt(h); 
% 
%         numeratore = drift_corrente * h + X_corrente - X_futuro_giu;
%         denominatore = 2 * sqrt(h);
% 
%         probabilita = max(0, min(1, numeratore / denominatore));
% 
%         jd_S(i+1, 1:i+1, 1:i+1)    = salto_giu + 1;
%         ju_S(i+1, 1:i+1, 1:i+1)    = salto_su + 1;
%         p_hat_S(i+1, 1:i+1, 1:i+1) = probabilita;
%     end
% 
%     
%     S_tree = NaN(N+1, N+1, N+1);
%     nodi_temporali = (0:N) * h;
% 
%     for i = 0:N
%         tempo_i = i * h;
%         integrale_I = interp1(nodi_temporali, I, tempo_i, 'linear', 'extrap');
% 
%         X_corrente    = X_tree(i+1, 1:i+1)';
%         zeta_corrente = zeta_tree(i+1, 1:i+1);
% 
%         esponente = sigma_S * (sqrt(1-ro^2) * X_corrente + ro * zeta_corrente) + ...
%                     integrale_I - (sigma_S^2 / 2) * tempo_i;
% 
%         S_matrice = S0 * exp(esponente);
% 
%         S_tree(i+1, 1:i+1, 1:i+1) = S_matrice;
%     end
% 
% end