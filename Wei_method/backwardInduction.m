%% VECTORAIZED FUNCTION

function V0 = backwardInduction(S_tree, x_tree, N, K, p_r, p_hat_S, kd_r,...
    ku_r, jd_S, ju_S, df_HW, rho, sigma_r, sigma_S, h, method)
% BACKWARDINDUCTION Prices an American FX Option using a Bivariate Tree.
% Handles both the Standard 3D approach ('Wei') and the 2D Copula approach ('Robust').
% BACKWARDINDUCTION Prices an American FX Option using a Bivariate Tree.
% Handles both the Standard 3D approach ('Wei') and the 2D Copula approach ('Robust').
%
% INPUTS:
%   S_tree  - FX rate spatial grid. 3D array (Wei) or 2D matrix (Robust).
%   x_tree  - HW stochastic rate spatial grid. 2D matrix.
%   N       - Number of time steps.
%   K       - Strike price of the option.
%   p_r     - Up-jump probabilities for the HW rate tree. 2D matrix.
%   p_hat_S - Marginal up-jump probabilities for the FX tree. 3D array or 2D matrix.
%   kd_r    - Down-jump indices for the HW rate tree. 2D matrix.
%   ku_r    - Up-jump indices for the HW rate tree. 2D matrix.
%   jd_S    - Down-jump indices for the FX tree. 3D array or 2D matrix.
%   ju_S    - Up-jump indices for the FX tree. 3D array or 2D matrix.
%   df_HW   - Local discount factors from the HW rate tree. 2D matrix.
%   rho     - Correlation coefficient between FX and domestic rate.
%   sigma_r - Volatility of the domestic interest rate.
%   sigma_S - Volatility of the FX spot rate.
%   h       - Time step size (dt).
%   method  - String specifying the algorithm: 'Wei' or 'Robust'.
%
% OUTPUT:
%   V0      - Scalar. The present value of the American option at time t=0.

    % Initialize the 3D tensor for option values
    V = NaN(N+1, N+1, N+1);
    
    % Detect if the input FX tree is 2D (Robust method) or 3D (Wei method)
    isRobustMethod = ismatrix(S_tree); 
    
    %  STEP 1: Terminal Payoff at Maturity (T) VECTORIZED 
    % Eliminate the two nested loops (j and k) by calculating the entire V(N+1) slice at once
    if isRobustMethod
        % S_end becomes a column vector (N+1 x 1)
        S_end = S_tree(N+1, 1:N+1)'; 
        % Compute the payoff for each j and replicate it N+1 times along k 
        % to fill the (N+1 x N+1) matrix
        V_end_slice = repmat(max(K - S_end, 0), 1, N+1);
    else 
        % Directly extract the entire 2D matrix at the end of the tree
        S_end = squeeze(S_tree(N+1, 1:N+1, 1:N+1));
        V_end_slice = max(K - S_end, 0);
    end
    % Insert the slice into the tensor (using reshape to ensure 1xNxN dimensions)
    V(N+1, 1:N+1, 1:N+1) = reshape(V_end_slice, [1, N+1, N+1]);
    
    %  STEP 2: Backward Induction Loop 
    % Keep only the time step loop "i", computing the (i x i) matrices instantaneously
    for i = N : -1 : 1
        
        %  DATA EXTRACTION FOR STEP i 
        % 'squeeze' removes singleton dimensions, returning (i x i) matrices
        ps = squeeze(p_hat_S(i, 1:i, 1:i)); 
        jd = squeeze(jd_S(i, 1:i, 1:i));    
        ju = squeeze(ju_S(i, 1:i, 1:i));    
        
        % Rate-related variables (index k) are row vectors (1 x i). 
        % We 'repmat' them to create (i x i) matrices needed to map the joint grid
        pr_mat = repmat(p_r(i, 1:i), i, 1);          
        kd_mat = repmat(kd_r(i, 1:i), i, 1);         
        ku_mat = repmat(ku_r(i, 1:i), i, 1);         
        df_mat = repmat(df_HW(i, 1:i), i, 1);        

        % Extract the current spot price
        if isRobustMethod
            S_cur = S_tree(i, 1:i)'; % Column vector (i x 1) -> maps to index 'j'
        else
            S_cur = squeeze(S_tree(i, 1:i, 1:i)); % Matrix (i x i)
        end
        
        %  STEP 3: Compute Joint Probabilities 
        switch method
            case 'Wei'
                % Equation 10
                p_UU = ps .* pr_mat;
                p_UD = ps .* (1 - pr_mat);
                p_DU = (1 - ps) .* pr_mat;
                p_DD = (1 - ps) .* (1 - pr_mat);
                
            case 'Robust'
                % Get necessary variables

                % FX Side
                S_next = S_tree(i+1, 1:i+1); % Row vector 1 x (i+1)
                S_u = S_next(ju); % Result: Matrix (i x i)
                S_d = S_next(jd); % Result: Matrix (i x i)
                
                % Stochastic rate Side 
                x_c = x_tree(i, 1:i);           % Row vector 1 x i
                x_next = x_tree(i+1, 1:i+1);    % Row vector 1 x (i+1)
                x_u = x_next(ku_r(i, 1:i));     % Row vector 1 x i
                x_d = x_next(kd_r(i, 1:i));     % Row vector 1 x i
                
                % Equation 22
                m_uu = (S_u - S_cur) .* (x_u - x_c);
                m_ud = (S_u - S_cur) .* (x_d - x_c);
                m_du = (S_d - S_cur) .* (x_u - x_c);
                m_dd = (S_d - S_cur) .* (x_d - x_c);
                
                % Left hand side of the 4th equality of equation 21
                C = rho * sigma_r * sigma_S * S_cur * h;
                % Analytical solution for the linear system denominator (can be close to 0)
                D = m_uu - m_ud - m_du + m_dd;
                %One liner system (21) solver
                q1 = (C - m_ud .* ps - m_du .* pr_mat - m_dd .* (1 - ps - pr_mat)) ./ D;
                
                % FAILSAFE: Handling D == 0 via Logical Indexing
                % Create a boolean mask (true/false) for points where the division failed
                zero_D = (D == 0);
                if any(zero_D, 'all')
                    % Overwrite q1 ONLY where D == 0, falling back to local independence
                    q1(zero_D) = ps(zero_D) .* pr_mat(zero_D); 
                end
                
                %  VECTORIZED POSITIVITY CHECK 
                % These bounds ensure all 4 derived joint probabilities remain >= 0
                p_UU_min = max(0, ps + pr_mat - 1);
                p_UU_max = min(ps, pr_mat);
                %In paper denoted q but for downstream simplicity here it's p
                p_UU = max(p_UU_min, min(p_UU_max, q1));
                p_UD = ps - p_UU;
                p_DU = pr_mat - p_UU;
                p_DD = 1 - ps - pr_mat + p_UU;
        end
        
        %  STEP 4: Discounting & Early Exercise 
        % Extract the option slice from the future time step i+1
        V_next = squeeze(V(i+1, 1:i+1, 1:i+1)); % Matrix (i+1 x i+1)
        
        % SUB2IND TRICK: Convert 2D jump indices [row, col] into a single linear index. 
        % This allows fetching values from V_next simultaneously for all (i x i) nodes, without loops!
        idx_UU = sub2ind([i+1, i+1], ju, ku_mat);
        idx_UD = sub2ind([i+1, i+1], ju, kd_mat);
        idx_DU = sub2ind([i+1, i+1], jd, ku_mat);
        idx_DD = sub2ind([i+1, i+1], jd, kd_mat);
        
        % Compute the expected future value
        V_UU_val = V_next(idx_UU) .* p_UU;
        V_UD_val = V_next(idx_UD) .* p_UD;
        V_DU_val = V_next(idx_DU) .* p_DU;
        V_DD_val = V_next(idx_DD) .* p_DD;
        
        % Discount using the local market discount factor df_HW (contains full r_t)
        continuationValue = df_mat .* (V_UU_val + V_UD_val + V_DU_val + V_DD_val);
        
        % Apply American early exercise boundary constraint
       
        V_curr = max(max(K - S_cur, 0), continuationValue);
        
        % Save results into the 3D tensor
        V(i, 1:i, 1:i) = reshape(V_curr, [1, i, i]);
    end
    
    % Final Price
    V0 = V(1, 1, 1);
end

% Legacy non vecto function
% function V0 = backwardInduction(S_tree, x_tree, N, K, p_r, p_hat_S,...
%  kd_r, ku_r, jd_S, ju_S, df_HW, rho, sigma_r, sigma_S, h, method)
% % BACKWARDINDUCTION Prices an American FX Option using a Bivariate Tree.
% % Handles both the Standard 3D approach ('Wei') and the 2D Copula approach ('Robust').
% 
%     % Initialize the 3D tensor for option values
%     V = NaN(N+1, N+1, N+1);
% 
%     % Detect if the input FX tree is 2D (Robust method) or 3D (Wei method)
%     isRobustMethod = ismatrix(S_tree); 
% 
%     %  STEP 1: Terminal Payoff at Maturity (T) 
%     for j = 1 : N+1
%         for k = 1: N+1
%             if isRobustMethod
%                 S_end = S_tree(N+1, j);
%             else 
%                 S_end = S_tree(N+1, j, k);
%             end
%             % Option intrinsic value at maturity
%             V(N+1, j, k) = max(K - S_end, 0);
%         end
%     end
% 
%     %  STEP 2: Backward Induction Loop 
%     for i = N : -1 : 1
%         for j = 1 : i
%             for k = 1 : i
% 
%                 % Extract marginal probabilities for the current node
%                 ps = p_hat_S(i, j, k);
%                 pr = p_r(i, k);
% 
%                 % Extract spatial jump indices
%                 jd = jd_S(i, j, k);
%                 ju = ju_S(i, j, k);
%                 kd = kd_r(i, k);
%                 ku = ku_r(i, k);
% 
%                 % Extract the current spot price
%                 if isRobustMethod
%                     S_cur = S_tree(i, j);
%                 else
%                     S_cur = S_tree(i, j, k);
%                 end
% 
%                 %  STEP 3: Compute Joint Probabilities 
%                 switch method
%                     case 'Wei'
%                         % Assumption of zero correlation: joint prob is the product of marginals
%                         p_UU = ps * pr;
%                         p_UD = ps * (1 - pr);
%                         p_DU = (1 - ps) * pr;
%                         p_DD = (1 - ps) * (1 - pr);
% 
%                     case 'Robust'
%                         % Extract future spot prices
%                         S_u = S_tree(i+1, ju);
%                         S_d = S_tree(i+1, jd);
% 
%                         % Extract pure stochastic rate process (x_t) for exact covariance matching
%                         x_c = x_tree(i, k);
%                         x_u = x_tree(i+1, ku);
%                         x_d = x_tree(i+1, kd);
% 
%                         % Cross-moments based purely on the stochastic component (isolating volatility)
%                         m_uu = (S_u - S_cur) * (x_u - x_c);
%                         m_ud = (S_u - S_cur) * (x_d - x_c);
%                         m_du = (S_d - S_cur) * (x_u - x_c);
%                         m_dd = (S_d - S_cur) * (x_d - x_c);
% 
%                         % Target covariance from input correlation
%                         C = rho * sigma_r * sigma_S * S_cur * h;
% 
%                         % Analytical solution for the linear system denominator
%                         D = m_uu - m_ud - m_du + m_dd;
% 
%                         % Compute unbounded (raw) joint probability
%                         if D ~= 0
%                             q1 = (C - m_ud*ps - m_du*pr - m_dd*(1 - ps - pr)) / D;
%                         % FAILSAFE: If D == 0, future nodes have collapsed (S_u == S_d or x_u == x_d),
%                         % risking a fatal division by zero. Since geometric correlation cannot be 
%                         % expressed in a collapsed grid, we fallback to local independence 
%                         % (zero correlation), where the joint probability is the product of the marginals.
%                         else
%                             q1 = ps * pr; 
%                         end
% 
%                         %  POSITIVITY CHECK: x
%                         % These bounds ensure all 4 derived joint probabilities remain >= 0
%                         % without destroying the marginal probabilities (which govern the risk-neutral drift).
%                         p_UU_min = max(0, ps + pr - 1);
%                         p_UU_max = min(ps, pr);
% 
%                         % Strict Clamping to physical correlation limits
%                         p_UU = max(p_UU_min, min(p_UU_max, q1));
% 
%                         % Cascade derivation of remaining probabilities
%                         % By clamping p_UU, these are mathematically guaranteed to be >= 0 
%                         % and to sum exactly to 1.0.
%                         p_UD = ps - p_UU;
%                         p_DU = pr - p_UU;
%                         p_DD = 1 - ps - pr + p_UU;
% 
%                 end
% 
%                 %  STEP 4: Discounting & Early Exercise 
%                 % Compute the expected future value of the option
%                 V_UU = V(i+1, ju, ku) * p_UU;
%                 V_UD = V(i+1, ju, kd) * p_UD;
%                 V_DU = V(i+1, jd, ku) * p_DU;
%                 V_DD = V(i+1, jd, kd) * p_DD;
% 
%                 % Discount using the local market discount factor df_HW (contains full r_t)
%                 continuationValue = df_HW(i, k) * (V_UU + V_UD + V_DU + V_DD);
% 
%                 % Apply American early exercise boundary constraint
%                 V(i, j, k) = max(max(K - S_cur, 0), continuationValue);
%             end
%         end
%     end
% 
%     % The root node contains the present value of the option
%     V0 = V(1, 1, 1);
% end
