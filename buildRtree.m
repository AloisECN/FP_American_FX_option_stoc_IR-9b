function [r_tree, zeta_tree, kd_r, ku_r, p_r] = buildRtree(N, h, kappa, sigma_r, zeta0)

    zeta_tree = zeros(N+1, N+1);
    for i = 0:N
        for k = 0:i
           zeta_tree(i+1, k+1) = zeta0 + (2*k - i) * sqrt(h);
        end
    end

    kd_r = NaN(N+1, N+1);
    ku_r = NaN(N+1, N+1);
    p_r  = NaN(N+1, N+1);

    for i = 0:N-1
        for k = 0:i
            drift = -kappa * zeta_tree(i+1, k+1);
            kd_r(i+1, k+1) = k + floor((drift * sqrt(h) + 1)/2);
            kd_r(i+1, k+1) = max(0, min(i, kd_r(i+1, k+1)));
            ku_r(i+1,k+1) = kd_r(i+1,k+1) + 1;

            num = drift * h + zeta_tree(i+1, k+1) - zeta_tree(i+2, kd_r(i+1, k+1)+1);
            den = 2*sqrt(h);
            p_r(i+1, k+1) = max(0, min(1, num/den));
        end

    end

    r_tree = sigma_r * zeta_tree;

end




%% VECTORIALIZED VERSION -> MORE EFFICIENT


% function [r_tree, kd_r, ku_r, p_r] = buildRtree_opt(N, h, kappa, sigma_r, zeta0)
% 
%     
%     i_vec = (0:N)'; 
%     k_vec = (0:N);  
% 
%     
%     kd_r = NaN(N+1, N+1);
%     ku_r = NaN(N+1, N+1);
%     p_r  = NaN(N+1, N+1);
% 
%     zeta_tree = zeta0 + (2 * k_vec - i_vec) * sqrt(h);
% 
%     valid_mask = k_vec <= i_vec;
%     zeta_tree(~valid_mask) = 0; 
% 
%     i_sub = i_vec(1:N); 
%     zeta_sub = zeta_tree(1:N, :);
% 
%     drift = -kappa * zeta_sub;
% 
%     kd_r_sub = k_vec + floor((drift * sqrt(h) + 1) / 2);
% 
%     kd_r_sub = max(0, min(i_sub, kd_r_sub));
%     ku_r_sub = kd_r_sub + 1;
% 
%     zeta_next_down = zeta0 + (2 * kd_r_sub - (i_sub + 1)) * sqrt(h);
% 
%     num = drift * h + zeta_sub - zeta_next_down;
%     den = 2 * sqrt(h); % Semplificazione matematica esatta del tuo denominatore
% 
%     p_r_sub = max(0, min(1, num ./ den));
% 
%     mask_sub = valid_mask(1:N, :);
%     kd_r_sub(~mask_sub) = NaN;
%     ku_r_sub(~mask_sub) = NaN;
%     p_r_sub(~mask_sub)  = NaN;
% 
%     kd_r(1:N, :) = kd_r_sub;
%     ku_r(1:N, :) = ku_r_sub;
%     p_r(1:N, :)  = p_r_sub;
% 
%     r_tree = sigma_r * zeta_tree;
% 
% end