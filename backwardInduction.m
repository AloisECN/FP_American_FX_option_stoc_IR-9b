% function V0 = backwardInduction(S_tree, N, K, p_r, p_hat_S, kd_r, jd_S, df_HW)
% 
% V = NaN(N+1, N+1, N+1);
% for j = 1 : N+1
%     for k = 1: N+1
%         V(N+1, j, k) = max(K - S_tree(N+1, j, k), 0);
%     end
% end
% 
% for i = N : -1 : 1
%     for j = 1 : i
%         for k = 1 : i
%             ps = p_hat_S(i, j, k);
%             pr = p_r(i, k);
% 
%             p_UU = ps * pr;
%             p_UD = ps * (1 - pr);
%             p_DU = (1 - ps) * pr;
%             p_DD = (1 - ps) * (1 - pr);
% 
%             jd = jd_S(i, j, k);
%             ju = 1 + jd;
%             kd = kd_r(i, k);
%             ku = 1 + kd;
% 
%             V_UU = V(i+1, ju, ku) * p_UU;
%             V_UD = V(i+1, ju, kd) * p_UD;
%             V_DU = V(i+1, jd, ku) * p_DU;
%             V_DD = V(i+1, jd, kd) * p_DD;
% 
%             df = df_HW(i, k);
%             continuationValue = df * (V_UU + V_UD + V_DU + V_DD);
% 
%             intrinsicValue = max(K - S_tree(i, j, k), 0);
% 
%             V(i, j, k) = max(intrinsicValue, continuationValue);
% 
%         end
%     end
% end
% 
% 
% V0 = V(1, 1, 1);
% end


function V0 = backwardInduction(S_tree, N, K, p_r, p_hat_S, kd_r, jd_S, df_HW)

    V = NaN(N+1, N+1, N+1);

    % 1. Vettorializzazione del Payoff Terminale (Step N+1)
    % Usiamo squeeze per convertire 1x(N+1)x(N+1) in una matrice (N+1)x(N+1)
    S_terminal = squeeze(S_tree(N+1, 1:N+1, 1:N+1));
    V(N+1, 1:N+1, 1:N+1) = max(K - S_terminal, 0);

    % 2. Ciclo di Backward Induction (Scorre solo il tempo all'indietro)
    for i = N : -1 : 1
        
        % --- ESTRAZIONE VARIABILI CORRENTI ---
        % Matrici (i x i) che dipendono da j (righe) e k (colonne)
        ps_mat = squeeze(p_hat_S(i, 1:i, 1:i)); 
        jd_mat = squeeze(jd_S(i, 1:i, 1:i));
        ju_mat = jd_mat + 1;
        
        % Vettori riga (1 x i) che dipendono SOLO da k (colonne)
        % Forziamo la forma a riga con reshape per l'espansione implicita
        pr_row = reshape(p_r(i, 1:i), 1, i);
        kd_row = reshape(kd_r(i, 1:i), 1, i);
        ku_row = kd_row + 1;
        df_row = reshape(df_HW(i, 1:i), 1, i);

        % --- CALCOLO PROBABILITA' CONGIUNTE ---
        % Moltiplicando una matrice (i x i) per un vettore riga (1 x i),
        % MATLAB applica il vettore ad ogni riga della matrice automaticamente.
        p_UU = ps_mat .* pr_row;
        p_UD = ps_mat .* (1 - pr_row);
        p_DU = (1 - ps_mat) .* pr_row;
        p_DD = (1 - ps_mat) .* (1 - pr_row);

        % --- RECUPERO VALORI FUTURI TRAMITE INDICIZZAZIONE LINEARE ---
        % Estraiamo l'intera matrice dei valori allo step i+1
        V_next = squeeze(V(i+1, 1:i+1, 1:i+1)); % Dimensione: (i+1) x (i+1)
        
        % Numero di righe della matrice futura (serve per calcolare l'indice lineare)
        num_rows = i + 1;
        
        % Formula Indice Lineare = (Colonna - 1) * num_righe + Riga
        % Riga dipende da j (ju_mat, jd_mat). Colonna dipende da k (ku_row, kd_row)
        % L'espansione implicita qui creerà istantaneamente matrici (i x i)
        idx_UU = (ku_row - 1) * num_rows + ju_mat;
        idx_UD = (kd_row - 1) * num_rows + ju_mat;
        idx_DU = (ku_row - 1) * num_rows + jd_mat;
        idx_DD = (kd_row - 1) * num_rows + jd_mat;

        % Estraiamo i valori e li moltiplichiamo per le probabilità in un colpo solo
        V_UU = V_next(idx_UU) .* p_UU;
        V_UD = V_next(idx_UD) .* p_UD;
        V_DU = V_next(idx_DU) .* p_DU;
        V_DD = V_next(idx_DD) .* p_DD;

        % --- CONTINUATION VALUE & INTRINSIC VALUE ---
        % Scontiamo i valori futuri
        continuationValue = df_row .* (V_UU + V_UD + V_DU + V_DD);

        % Valore di esercizio anticipato
        S_curr = squeeze(S_tree(i, 1:i, 1:i));
        intrinsicValue = max(K - S_curr, 0);

        % Assegnazione del massimo (American Option logic) alla slice 3D
        V(i, 1:i, 1:i) = max(intrinsicValue, continuationValue);

    end

    % Il valore al tempo 0
    V0 = V(1, 1, 1);
end