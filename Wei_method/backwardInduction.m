function V0 = backwardInduction(S_tree, x_tree, N, K, p_r, p_hat_S, kd_r, ku_r, jd_S, ju_S, df_HW, rho, sigma_r, sigma_S, h, method)
% BACKWARDINDUCTION Prices an American FX Option using a Bivariate Tree.
% Handles both the Standard 3D approach ('Wei') and the 2D Copula approach ('Robust').

    % Initialize the 3D tensor for option values
    V = NaN(N+1, N+1, N+1);
    
    % Detect if the input FX tree is 2D (Robust method) or 3D (Wei method)
    isRobustMethod = ismatrix(S_tree); 
    
    % --- STEP 1: Terminal Payoff at Maturity (T) ---
    for j = 1 : N+1
        for k = 1: N+1
            if isRobustMethod
                S_end = S_tree(N+1, j);
            else 
                S_end = S_tree(N+1, j, k);
            end
            % Option intrinsic value at maturity
            V(N+1, j, k) = max(K - S_end, 0);
        end
    end
    
    % --- STEP 2: Backward Induction Loop ---
    for i = N : -1 : 1
        for j = 1 : i
            for k = 1 : i
                
                % Extract marginal probabilities for the current node
                ps = p_hat_S(i, j, k);
                pr = p_r(i, k);
                
                % Extract spatial jump indices
                jd = jd_S(i, j, k);
                ju = ju_S(i, j, k);
                kd = kd_r(i, k);
                ku = ku_r(i, k);
                
                % Extract the current spot price
                if isRobustMethod
                    S_cur = S_tree(i, j);
                else
                    S_cur = S_tree(i, j, k);
                end
                
                % --- STEP 3: Compute Joint Probabilities ---
                switch method
                    case 'Wei'
                        % Assumption of zero correlation: joint prob is the product of marginals
                        p_UU = ps * pr;
                        p_UD = ps * (1 - pr);
                        p_DU = (1 - ps) * pr;
                        p_DD = (1 - ps) * (1 - pr);
                        
                    case 'Robust'
                        % Extract future spot prices
                        S_u = S_tree(i+1, ju);
                        S_d = S_tree(i+1, jd);
                        
                        % Extract pure stochastic rate process (x_t) for exact covariance matching
                        x_c = x_tree(i, k);
                        x_u = x_tree(i+1, ku);
                        x_d = x_tree(i+1, kd);
                        
                        % Cross-moments based purely on the stochastic component (isolating volatility)
                        m_ud = (S_u - S_cur) * (x_d - x_c);
                        m_du = (S_d - S_cur) * (x_u - x_c);
                        m_dd = (S_d - S_cur) * (x_d - x_c);
                        
                        % Target covariance from input correlation
                        C = rho * sigma_r * sigma_S * S_cur * h;
                        
                        % Analytical solution for the linear system denominator
                        D = (S_u - S_d) * (x_u - x_d);
                        
                        % Compute unbounded (raw) joint probability
                        if D ~= 0
                            q1_raw = (C - m_ud*ps - m_du*pr - m_dd*(1 - ps - pr)) / D;
                        else
                            q1_raw = ps * pr; % Failsafe if nodes collapse
                        end
                        
                        % --- POSITIVITY CHECK: x---
                        % These bounds ensure all 4 derived joint probabilities remain >= 0
                        % without destroying the marginal probabilities (which govern the risk-neutral drift).
                        p_UU_min = max(0, ps + pr - 1);
                        p_UU_max = min(ps, pr);
                        
                        % Strict Clamping to physical correlation limits
                        p_UU = max(p_UU_min, min(p_UU_max, q1_raw));
                        
                        % Cascade derivation of remaining probabilities
                        % By clamping p_UU, these are mathematically guaranteed to be >= 0 
                        % and to sum exactly to 1.0.
                        p_UD = ps - p_UU;
                        p_DU = pr - p_UU;
                        p_DD = 1 - ps - pr + p_UU;
                end
                
                % --- STEP 4: Discounting & Early Exercise ---
                % Compute the expected future value of the option
                V_UU = V(i+1, ju, ku) * p_UU;
                V_UD = V(i+1, ju, kd) * p_UD;
                V_DU = V(i+1, jd, ku) * p_DU;
                V_DD = V(i+1, jd, kd) * p_DD;
                
                % Discount using the local market discount factor df_HW (contains full r_t)
                continuationValue = df_HW(i, k) * (V_UU + V_UD + V_DU + V_DD);
                
                % Apply American early exercise boundary constraint
                V(i, j, k) = max(max(K - S_cur, 0), continuationValue);
            end
        end
    end
    
    % The root node contains the present value of the option
    V0 = V(1, 1, 1);
end



% function V0 = backwardInduction(S_tree, x_tree, N, K, p_r, p_hat_S, kd_r, ku_r, jd_S, ju_S, df_HW, rho, sigma_r, sigma_S, h, method)
% 
%     V = NaN(N+1, N+1, N+1);
%     isRobustMethod = ismatrix(S_tree); 
% 
%     if isRobustMethod
%         S_end = S_tree(N+1, 1:N+1); 
%         payoff = max(K - S_end, 0); 
% 
%         payoff_3D = repmat(reshape(payoff, 1, N+1, 1), 1, 1, N+1);
%         V(N+1, :, :) = payoff_3D; 
%     else 
%         S_end = S_tree(N+1, 1:N+1, 1:N+1);
%         V(N+1, :, :) = max(K - S_end, 0);
%     end
% 
%     for i = N : -1 : 1
% 
%         ps = squeeze(p_hat_S(i, 1:i, 1:i)); 
%         pr = p_r(i, 1:i);                   
% 
%         jd = squeeze(jd_S(i, 1:i, 1:i));    
%         ju = squeeze(ju_S(i, 1:i, 1:i));    
% 
%         kd = kd_r(i, 1:i);                  
%         ku = ku_r(i, 1:i);                       
% 
%         if isRobustMethod
%             S_cur = S_tree(i, 1:i)';        
%         else
%             S_cur = squeeze(S_tree(i, 1:i, 1:i)); 
%         end
% 
%         switch method
%             case 'Wei'
%                 p_UU = ps .* pr;
%                 p_UD = ps .* (1 - pr);
%                 p_DU = (1 - ps) .* pr;
%                 p_DD = (1 - ps) .* (1 - pr);
% 
%             case 'Robust'
%                 S_next = S_tree(i+1, :);       
%                 S_u = S_next(ju);              
%                 S_d = S_next(jd);              
% 
%                 % --- CORREZIONE VETTORIALIZZATA: Tassi Pieni dagli Sconti ---
%                 r_c = -log(df_HW(i, 1:i)) / h;       
% 
%                 df_next = df_HW(i+1, 1:i+1);
%                 r_next_full = -log(df_next) / h;
%                 r_u = r_next_full(ku);      
%                 r_d = r_next_full(kd);      
% 
%                 % Il resto del calcolo tensoriale rimane identico e perfettamente funzionante
%                 m_uu = (S_u - S_cur) .* (r_u - r_c);
%                 m_ud = (S_u - S_cur) .* (r_d - r_c);
%                 m_du = (S_d - S_cur) .* (r_u - r_c);
%                 m_dd = (S_d - S_cur) .* (r_d - r_c);
% 
%                 C = rho * sigma_r * sigma_S * S_cur * h; 
% 
%                 Denom = m_uu - m_ud - m_du + m_dd;
%                 Denom(Denom == 0) = eps; 
% 
%                 q1 = (C - m_ud .* ps - m_du .* pr - m_dd .* (1 - ps - pr)) ./ Denom;
%                 q2 = ps - q1;
%                 q3 = pr - q1;
%                 q4 = 1 - ps - pr + q1;
% 
%                 p_UU_raw = max(0, q1);
%                 p_UD_raw = max(0, q2);
%                 p_DU_raw = max(0, q3);
%                 p_DD_raw = max(0, q4);
% 
%                 sum_p = p_UU_raw + p_UD_raw + p_DU_raw + p_DD_raw;
%                 mask = sum_p > 0;
% 
%                 p_UU = repmat(0.25, i, i); p_UD = p_UU; p_DU = p_UU; p_DD = p_UU;
% 
%                 p_UU(mask) = p_UU_raw(mask) ./ sum_p(mask);
%                 p_UD(mask) = p_UD_raw(mask) ./ sum_p(mask);
%                 p_DU(mask) = p_DU_raw(mask) ./ sum_p(mask);
%                 p_DD(mask) = p_DD_raw(mask) ./ sum_p(mask);
%         end
% 
%         V_next = squeeze(V(i+1, 1:i+1, 1:i+1));
% 
%         KU = repmat(ku, i, 1);
%         KD = repmat(kd, i, 1);
% 
%         idx_uu = sub2ind([i+1, i+1], ju, KU);
%         idx_ud = sub2ind([i+1, i+1], ju, KD);
%         idx_du = sub2ind([i+1, i+1], jd, KU);
%         idx_dd = sub2ind([i+1, i+1], jd, KD);
% 
%         V_UU = V_next(idx_uu) .* p_UU;
%         V_UD = V_next(idx_ud) .* p_UD;
%         V_DU = V_next(idx_du) .* p_DU;
%         V_DD = V_next(idx_dd) .* p_DD;
% 
%         df = df_HW(i, 1:i); 
% 
%         continuationValue = df .* (V_UU + V_UD + V_DU + V_DD); 
% 
%         V_current = max(max(K - S_cur, 0), continuationValue);
% 
%         V(i, 1:i, 1:i) = reshape(V_current, 1, i, i);
% 
%     end
% 
%     V0 = V(1, 1, 1);
% end