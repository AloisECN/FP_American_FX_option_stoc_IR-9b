%% Vectorialized function

function [S_tree_rob, jd_S_rob, ju_S_rob, p_hat_S_rob] = buildStreeRobust(N, h, sigma_S, U0, df_HW)
% BUILDSTREEROBUST Constructs the static FX spatial grid and computes 
% transition probabilities for the Robust Bivariate Tree (Appolloni method).
%
% INPUTS:
% N, h, sigma_S, U0: Standard tree parameters (steps, dt, vol, initial state)
% df_HW: The Hull-White local discount factor matrix (used to extract full rates)

    %  1. Pure Volatility Grid VECTORIZED 
    % Unlike the standard Wei method, this grid is "rigid" and decoupled from 
    % the interest rate. It expands symmetrically based purely on the FX 
    % continuous volatility (sigma_S), avoiding distortion.
    
    i_grid = (0:N)';   
    j_grid = 0:N;      
    %Equation 16
    U_val = U0 + (2*j_grid - i_grid) * sqrt(h); 
    
    S_tree_rob = zeros(N+1, N+1);
    valid_mask = (j_grid <= i_grid);
    %Equation 17
    S_tree_rob(valid_mask) = exp(sigma_S * U_val(valid_mask));

    jd_S_rob = NaN(N+1, N+1, N+1);
    ju_S_rob = NaN(N+1, N+1, N+1);
    p_hat_S_rob = NaN(N+1, N+1, N+1);

    %Semi-vectorised we still loop through the i bu do the N+1xN+1 matrix 
    for i = 1 : N
        %get constants for future computations
        r_t = -log(df_HW(i, 1:i)) / h;
        S_cur = S_tree_rob(i, 1:i)';
        S_next = S_tree_rob(i+1, 1:i+1);
        %Infer it from our model (! diferent from the paper !)
        drift = S_cur .* r_t;
        threshold = S_cur + drift * h;

        % Since our grid is rigid, the 'threshold' falls between two predefined nodes.
        % Iterative 'find' is extremely slow on matrices. We use 'discretize' (binning)
        % to simultaneously map all (i x i) thresholds into their corresponding intervals.
        edges = [-Inf, S_next, Inf];
        jd_mat = discretize(threshold, edges) - 1; 
        %equation 18 & 19 vectorized

        jd_mat = max(1, min(i+1, jd_mat)); 
        ju_mat = min(i+1, jd_mat + 1);
        % Line below hadles ju == jd edge case
        ju_mat(threshold == S_next(jd_mat)) = jd_mat(threshold == S_next(jd_mat));
        %Line below handles if the treshold is too low for the available nodes
        ju_mat(threshold < S_next(1)) = 1;    
         %Line below analgously handle the too high problm
        ju_mat(threshold > S_next(end)) = i+1;   

        % Fetch the physical FX values corresponding to the jump indices
        S_jd = S_next(jd_mat);
        S_ju = S_next(ju_mat);
        num = threshold - S_jd;
        den = S_ju - S_jd;
        
        % Preallocate the probability matrix with the dummy 0.5 probability.
        % This automatically handles the collapse failsafe where ju == jd, 
        % preventing division by zero
        p_hat_mat = repmat(0.5, i, i);
        
        % Create a boolean mask where ceiling and floor are distinct nodes
        valid_jump = (ju_mat > jd_mat);
        
        p_hat_mat(valid_jump) = num(valid_jump) ./ den(valid_jump);
        
        % Trasfrom the prob to acceptable ones
        p_hat_mat = max(0, min(1, p_hat_mat));
        
        jd_S_rob(i, 1:i, 1:i) = reshape(jd_mat, [1, i, i]);
        ju_S_rob(i, 1:i, 1:i) = reshape(ju_mat, [1, i, i]);
        p_hat_S_rob(i, 1:i, 1:i) = reshape(p_hat_mat, [1, i, i]);

    end
end

%Legacy code of the non vectorised equivalent function

% function [S_tree_rob, jd_S_rob, ju_S_rob, p_hat_S_rob] = buildStreeRobust(N, h, sigma_S, U0, df_HW)
% % BUILDSTREEROBUST Constructs the static FX spatial grid and computes 
% % transition probabilities for the Robust Bivariate Tree (Appolloni method).
% %
% % INPUTS:
% % N, h, sigma_S, U0: Standard tree parameters (steps, dt, vol, initial state)
% % df_HW: The Hull-White local discount factor matrix (used to extract full rates)
% 
%     % Initialize the rigid spatial grid
%     S_tree_rob = zeros(N+1, N+1);
% 
%     %  Pure Volatility Grid 
%     % Unlike the standard Wei method, this grid is "rigid" and decoupled from 
%     % the interest rate. It expands symmetrically based purely on the FX 
%     % continuous volatility (sigma_S), avoiding distortion.
%     for i = 0:N
%         for j = 0:i
%             U_val = U0 + (2*j - i) * sqrt(h);
%             S_tree_rob(i+1, j+1) = exp(sigma_S * U_val);
%         end
%     end
% 
%     % Initialize output tensors
%     jd_S_rob = NaN(N+1, N+1, N+1);
%     ju_S_rob = NaN(N+1, N+1, N+1);
%     p_hat_S_rob = NaN(N+1, N+1, N+1);
% 
%     %  Main Loop: Computing Transition Dynamics 
%     for i = 1 : N
%         for j = 1 : i
%             for k = 1 : i
% 
%                 %  Full Market Rate Extraction 
%                 % To maintain exact risk-neutral pricing and align with the Wei method, 
%                 % the drift must be driven by the FULL market rate (r_t), not just the 
%                 % stochastic pure process (x_t). We extract this continuous rate 
%                 % directly from the local discount factor matrix.
%                 df_local = df_HW(i, k);
%                 r_t = -log(df_local) / h;
% 
%                 % Calculate the exact risk-neutral expected future value (continuous)
%                 drift = r_t * S_tree_rob(i, j);
%                 threshold = S_tree_rob(i, j) + drift * h;
% 
%                 %  Continuous vs Discrete Mapping 
%                 % Since our grid is rigid, the 'threshold' falls between two predefined 
%                 % nodes. We search the t+1 grid to find the immediate lower floor (jd) 
%                 % and immediate upper ceiling (ju) nodes.
%                 jd = find(S_tree_rob(i+1, 1:i+1) <= threshold, 1, 'last');
%                 ju = find(S_tree_rob(i+1, 1:i+1) >= threshold, 1, 'first');
% 
%                 %  Boundary Clamping (Tail Truncation) 
%                 % If an extreme interest rate shock pushes the expected FX value 
%                 % completely outside the physical boundaries of our rigid grid,
%                 % the 'find' function returns empty. We clamp the target to the 
%                 % outermost nodes to prevent indexing crashes and contain the probability mass.
%                 if isempty(jd), jd = 1; end
%                 if isempty(ju), ju = i+1; end
% 
%                 jd_S_rob(i, j, k) = jd;
%                 ju_S_rob(i, j, k) = ju;
% 
%                 %  The Collapse Failsafe (ju = jd) 
%                 if ju > jd
%                     % Standard linear interpolation to center the expected value
%                     num = threshold - S_tree_rob(i+1, jd);
%                     den = S_tree_rob(i+1, ju) - S_tree_rob(i+1, jd);
%                     p_hat_S_rob(i, j, k) = max(0, min(1, num/den));
%                 else
%                     % If the grid is too sparse, or if Boundary Clamping forced ju == jd, 
%                     % the ceiling and floor collapse into the exact same physical node.
%                     % We assign a dummy 0.5 probability to prevent a division by zero 
%                     % (den = 0). This preserves mathematical accuracy since E[V] = V_down 
%                     % regardless of the probability value assigned here.
%                     p_hat_S_rob(i, j, k) = 0.5;
%                 end
% 
%             end
%         end
%     end
% end
