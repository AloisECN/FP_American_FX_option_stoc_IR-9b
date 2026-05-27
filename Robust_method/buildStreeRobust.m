% function [S_tree_rob, jd_S_rob, ju_S_rob, p_hat_S_rob] = buildStreeRobust(N, h, sigma_S, U0, df_HW)
%     S_tree_rob = zeros(N+1, N+1);
% 
%     for i = 0:N
%         for j = 0:i
%             U_val = U0 + (2*j - i) * sqrt(h);
%             S_tree_rob(i+1, j+1) = exp(sigma_S * U_val);
%         end
%     end
% 
%     jd_S_rob = NaN(N+1, N+1, N+1);
%     ju_S_rob = NaN(N+1, N+1, N+1);
%     p_hat_S_rob = NaN(N+1, N+1, N+1);
% 
%     for i = 1 : N
%         for j = 1 : i
%             for k = 1 : i
%                 df_local = df_HW(i, k);
%                 r_t = -log(df_local) / h;
% 
%                 drift = r_t * S_tree_rob(i, j);
%                 threshold = S_tree_rob(i, j) + drift * h;
% 
%                 jd = find(S_tree_rob(i+1, 1:i+1) <= threshold, 1, 'last');
%                 ju = find(S_tree_rob(i+1, 1:i+1) >= threshold, 1, 'first');
% 
%                 if isempty(jd), jd = 1; end
%                 if isempty(ju), ju = i+1; end
% 
%                 jd_S_rob(i, j, k) = jd;
%                 ju_S_rob(i, j, k) = ju;
% 
%                 if ju > jd
%                     num = threshold - S_tree_rob(i+1, jd);
%                     den = S_tree_rob(i+1, ju) - S_tree_rob(i+1, jd);
%                     p_hat_S_rob(i, j, k) = max(0, min(1, num/den));
%                 else
%                     p_hat_S_rob(i, j, k) = 0.5;
%                 end
%             end
%         end
%     end
% end


function [S_tree_rob, jd_S_rob, ju_S_rob, p_hat_S_rob] = buildStreeRobust(N, h, sigma_S, U0, df_HW)
    
    i_grid = (0:N)';    
    j_grid = 0:N;       
    
    U_val = U0 + (2*j_grid - i_grid) * sqrt(h);
    
    S_tree_rob = zeros(N+1, N+1);
    mask_tree = (j_grid <= i_grid); 
    S_tree_rob(mask_tree) = exp(sigma_S * U_val(mask_tree));
    
    jd_S_rob = NaN(N+1, N+1, N+1);
    ju_S_rob = NaN(N+1, N+1, N+1);
    p_hat_S_rob = NaN(N+1, N+1, N+1);
    
    for i = 1 : N
        
        S_curr = S_tree_rob(i, 1:i)';      
        df_curr = df_HW(i, 1:i);           
        
        r_t_curr = -log(df_curr) / h;     
        drift = S_curr .* r_t_curr;        
        threshold = S_curr + drift * h;    
        
        S_next = S_tree_rob(i+1, 1:i+1);   
        
        edges = [-Inf, S_next, Inf];
        jd_mat = discretize(threshold, edges) - 1;
        
        jd_mat = max(1, min(i+1, jd_mat));
        
        ju_mat = min(i+1, jd_mat + 1);
        
        ju_mat(threshold == S_next(jd_mat)) = jd_mat(threshold == S_next(jd_mat)); 
        ju_mat(threshold < S_next(1)) = 1;      
        ju_mat(threshold > S_next(end)) = i+1;  
        
        S_jd = S_next(jd_mat); 
        S_ju = S_next(ju_mat); 
        
        num = threshold - S_jd;
        den = S_ju - S_jd;
        
        p_hat_mat = zeros(i, i);
        mask_diff = (ju_mat > jd_mat);
        
        p_hat_mat(mask_diff) = num(mask_diff) ./ den(mask_diff);
        p_hat_mat(~mask_diff) = 0.5;
        
        p_hat_mat = max(0, min(1, p_hat_mat));
        
        jd_S_rob(i, 1:i, 1:i) = reshape(jd_mat, [1, i, i]);
        ju_S_rob(i, 1:i, 1:i) = reshape(ju_mat, [1, i, i]);
        p_hat_S_rob(i, 1:i, 1:i) = reshape(p_hat_mat, [1, i, i]);
        
    end
end
