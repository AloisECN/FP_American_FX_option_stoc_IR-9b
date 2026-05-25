% function [x_tree, zeta_tree, kd_r, ku_r, p_r] = buildXtree(N, h, kappa, sigma_r, zeta0)
% 
%     zeta_tree = zeros(N+1, N+1);
%     for i = 1 : N+1
%         for k = 1 : i
%            zeta_tree(i, k) = zeta0 + (2*k - i - 1) * sqrt(h);
%         end
%     end
% 
%     kd_r = NaN(N+1, N+1);
%     ku_r = NaN(N+1, N+1);
%     p_r  = NaN(N+1, N+1);
% 
%     for i = 1 : N
%         for k = 1 : i
%             drift = -kappa * zeta_tree(i, k);
%             kd_r(i, k) = k + floor((drift * sqrt(h) + 1)/2);
% 
%             % 1 is the lowest index of the tree
%             kd_r(i, k) = max(1, min(i, kd_r(i, k)));
%             ku_r(i,k) = kd_r(i,k) + 1;
% 
%             num = drift * h + zeta_tree(i, k) - zeta_tree(i+1, kd_r(i, k));
%             den = 2 * sqrt(h);
%             p_r(i, k) = max(0, min(1, num/den));
%         end
% 
%     end
% 
%     x_tree = sigma_r * zeta_tree;
% 
% end
% 



function [x_tree, zeta_tree, kd_r, ku_r, p_r] = buildXtree(N, h, kappa, sigma_r, zeta0)

    % 1. Creiamo matrici di indici per i e k
    % I rappresenta le righe, K rappresenta le colonne
    [K, I] = meshgrid(1:N+1, 1:N+1);

    % 2. Vettorializzazione di zeta_tree
    zeta_tree = zeros(N+1, N+1);
    % Maschera per il triangolo inferiore (k <= i)
    mask_zeta = K <= I; 
    zeta_tree(mask_zeta) = zeta0 + (2*K(mask_zeta) - I(mask_zeta) - 1) * sqrt(h);

    % 3. Inizializzazione delle matrici di output
    kd_r = NaN(N+1, N+1);
    ku_r = NaN(N+1, N+1);
    p_r  = NaN(N+1, N+1);

    % 4. Vettorializzazione del blocco principale
    % Maschera equivalente a: for i = 1:N, for k = 1:i
    mask_loop = (I <= N) & (K <= I);

    % Estraiamo i vettori degli elementi validi per elaborare tutto in un colpo solo
    i_vec = I(mask_loop);
    k_vec = K(mask_loop);
    zeta_vec = zeta_tree(mask_loop);

    % Calcolo di drift e kd_r
    drift = -kappa * zeta_vec;
    kd_vec = k_vec + floor((drift * sqrt(h) + 1) / 2);
    
    % Applicazione dei limiti max e min
    kd_vec = max(1, min(i_vec, kd_vec));
    ku_vec = kd_vec + 1;

    % Calcolo per zeta_tree(i+1, kd_r(i, k)) usando l'indicizzazione lineare
    % Formula dell'indice lineare per una matrice di dimensione (N+1)x(N+1):
    % indice = (colonna - 1) * num_righe + riga
    lin_idx_next = (kd_vec - 1) * (N + 1) + (i_vec + 1);
    zeta_next = zeta_tree(lin_idx_next);

    % Calcolo delle probabilità
    num = drift * h + zeta_vec - zeta_next;
    den = 2 * sqrt(h);
    p_vec = max(0, min(1, num / den));

    % 5. Ricostruzione delle matrici finali riassegnando i vettori alle posizioni corrette
    kd_r(mask_loop) = kd_vec;
    ku_r(mask_loop) = ku_vec;
    p_r(mask_loop)  = p_vec;

    % 6. Calcolo finale x_tree (già vettorializzato)
    x_tree = sigma_r * zeta_tree;

end