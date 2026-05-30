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
    
    U_val = U0 + (2*j_grid - i_grid) * sqrt(h);
    
    S_tree_rob = zeros(N+1, N+1);
    valid_mask = (j_grid <= i_grid);
    S_tree_rob(valid_mask) = exp(sigma_S * U_val(valid_mask));

    jd_S_rob = NaN(N+1, N+1, N+1);
    ju_S_rob = NaN(N+1, N+1, N+1);
    p_hat_S_rob = NaN(N+1, N+1, N+1);

    for i = 1 : N
        
        % Extract the continuous rate directly from the local discount factor matrix.
        r_t = -log(df_HW(i, 1:i)) / h;

        % Calculate the exact risk-neutral expected future value (continuous).
        S_cur = S_tree_rob(i, 1:i)';
        
        % Implicit Expansion: Multiplying a Column (i x 1) by a Row (1 x i) 
        % instantly generates the full (i x i) matrices for drift and threshold.
        drift = S_cur .* r_t;
        threshold = S_cur + drift * h;

        % Extract the future physical nodes (1D array)
        S_next = S_tree_rob(i+1, 1:i+1);

        % Since our grid is rigid, the 'threshold' falls between two predefined nodes.
        % Iterative 'find' is extremely slow on matrices. We use 'discretize' (binning)
        % to simultaneously map all (i x i) thresholds into their corresponding intervals.
        edges = [-Inf, S_next, Inf];
        jd_mat = discretize(threshold, edges) - 1;

        %  Boundary Clamping (Tail Truncation) 
        % If an extreme interest rate shock pushes the expected FX value outside 
        % the physical boundaries, we clamp the target to the outermost nodes.
        jd_mat = max(1, min(i+1, jd_mat));
        ju_mat = min(i+1, jd_mat + 1);
        
        % Exact edge-case corrections (replicating the original 'find' behavior)
        ju_mat(threshold == S_next(jd_mat)) = jd_mat(threshold == S_next(jd_mat)); 
        ju_mat(threshold < S_next(1)) = 1;       
        ju_mat(threshold > S_next(end)) = i+1;   

        %  The Collapse Failsafe (ju = jd) 
        % Fetch the physical FX values corresponding to the jump indices
        S_jd = S_next(jd_mat);
        S_ju = S_next(ju_mat);
        
        num = threshold - S_jd;
        den = S_ju - S_jd;
        
        % Preallocate the probability matrix with the dummy 0.5 probability.
        % This automatically handles the collapse failsafe where ju == jd, 
        % preventing division by zero while preserving E[V] = V_down.
        p_hat_mat = repmat(0.5, i, i);
        
        % Create a boolean mask where ceiling and floor are distinct nodes
        valid_jump = (ju_mat > jd_mat);
        
        p_hat_mat(valid_jump) = num(valid_jump) ./ den(valid_jump);
        
        % Clamp probabilities to [0, 1]
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
