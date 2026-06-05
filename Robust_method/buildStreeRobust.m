%% Vectorialized function

function [S_tree_rob, jd_S_rob, ju_S_rob, p_hat_S_rob] = buildStreeRobust(N, h, sigma_S, U0, x_tree, I, yearFrac, shift_amount)
% BUILDSTREEROBUST Constructs the static FX spatial grid and computes 
% transition probabilities for the Robust Bivariate Tree (Appolloni method).
%
% INPUTS:
% N, h, sigma_S, U0: Standard tree parameters and initial states
% x_tree: 2D matrix of the stochastic component of the domestic rate
% I: Vector of the deterministic market forward curve integration components
% yearFrac: Time structure vector associated with the market curve I
% shift_amount: Optional parallel shift applied to the cost-of-carry (default 0)

    if nargin < 8
        shift_amount = 0;
    end

    %  1. Pure Volatility Grid VECTORIZED 
    % Unlike the standard Wei method, this grid is "rigid" and decoupled from 
    % the interest rate. It expands symmetrically based purely on the FX 
    % continuous volatility (sigma_S), avoiding distortion.
    
    i_grid = (0:N)';   
    j_grid = 0:N;      
    % Equation 16
    U_val = U0 + (2*j_grid - i_grid) * sqrt(h); 
    
    S_tree_rob = zeros(N+1, N+1);
    valid_mask = (j_grid <= i_grid);
    % Equation 17
    S_tree_rob(valid_mask) = exp(sigma_S * U_val(valid_mask));

    jd_S_rob = NaN(N+1, N+1, N+1);
    ju_S_rob = NaN(N+1, N+1, N+1);
    p_hat_S_rob = NaN(N+1, N+1, N+1);

    % Pre-interpolate the market integral I(t) to match lattice time-steps
    t_nodes = (0:N) * h;
    I_ti = interp1(yearFrac, I, t_nodes, 'linear', 'extrap');

    % Semi-vectorised: we loop through 'i' but compute the full (i x i) matrix simultaneously
    for i = 1 : N
        
        % Extract current and future FX spot nodes
        S_cur = S_tree_rob(i, 1:i)'; % Column vector (i x 1)
        S_next = S_tree_rob(i+1, 1:i+1);
        
        % Extract the stochastic rate directly from x_tree
        x_slice = x_tree(i, 1:i); % Row vector (1 x i)
        
        % Derive the local deterministic cost-of-carry from the integral
        c_local = (I_ti(i+1) - I_ti(i)) / h;
        
        % Compute the true FX risk-neutral drift: (r_t - d_t) = x_t + c_t + shift
        drift_rate = x_slice + c_local + shift_amount;
        
        % Implicit expansion: (i x 1) .* (1 x i) generates an (i x i) drift matrix
        drift = S_cur .* drift_rate;
        threshold = S_cur + drift * h;

        % Since our grid is rigid, the 'threshold' falls between two predefined nodes.
        % We use 'discretize' (binning) to map all thresholds into corresponding intervals.
        edges = [-Inf, S_next, Inf];
        jd_mat = discretize(threshold, edges) - 1; 
        
        % Equation 18 & 19 vectorized boundary enforcement
        jd_mat = max(1, min(i+1, jd_mat)); 
        ju_mat = min(i+1, jd_mat + 1);
        
        % Edge cases handling
        ju_mat(threshold == S_next(jd_mat)) = jd_mat(threshold == S_next(jd_mat));
        ju_mat(threshold < S_next(1)) = 1;    
        ju_mat(threshold > S_next(end)) = i+1;   

        % Fetch the physical FX values corresponding to the jump indices
        S_jd = S_next(jd_mat);
        S_ju = S_next(ju_mat);
        num = threshold - S_jd;
        den = S_ju - S_jd;
        
        % Preallocate the probability matrix with the 0.5 failsafe
        p_hat_mat = repmat(0.5, i, i);
        
        % Create a boolean mask where ceiling and floor are distinct nodes
        valid_jump = (ju_mat > jd_mat);
        
        % Compute exact probability using element-wise division for the matrix
        p_hat_mat(valid_jump) = num(valid_jump) ./ den(valid_jump);
        
        % Transform the prob to acceptable ones [0, 1]
        p_hat_mat = max(0, min(1, p_hat_mat));
        
        % Reshape and store in the 3D tensors
        jd_S_rob(i, 1:i, 1:i) = reshape(jd_mat, [1, i, i]);
        ju_S_rob(i, 1:i, 1:i) = reshape(ju_mat, [1, i, i]);
        p_hat_S_rob(i, 1:i, 1:i) = reshape(p_hat_mat, [1, i, i]);
    end
end

%%

% function [S_tree_rob, jd_S_rob, ju_S_rob, p_hat_S_rob] = buildStreeRobust_loops(N, h, sigma_S, U0, x_tree, I, yearFrac, shift_amount)
% % BUILDSTREEROBUST_LOOPS Constructs the static FX spatial grid and computes 
% % transition probabilities for the Robust Bivariate Tree (Appolloni method).
% % (NON-VECTORIZED LEGACY EQUIVALENT)
% %
% % INPUTS:
% % N, h, sigma_S, U0: Standard tree parameters and initial states
% % x_tree: 2D matrix of the stochastic component of the domestic rate
% % I: Vector of the deterministic market forward curve integration components
% % yearFrac: Time structure vector associated with the market curve I
% % shift_amount: Optional parallel shift applied to the cost-of-carry (default 0)
% 
%     if nargin < 8
%         shift_amount = 0;
%     end
% 
%     % Initialize the rigid spatial grid
%     S_tree_rob = zeros(N+1, N+1);
% 
%     % 1. Pure Volatility Grid (Loops)
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
%     % Pre-interpolate the market integral I(t) to match lattice time-steps
%     % (Done outside the loops to avoid redundant computations)
%     t_nodes = (0:N) * h;
%     I_ti = interp1(yearFrac, I, t_nodes, 'linear', 'extrap');
% 
%     % 2. Main Loop: Computing Transition Dynamics 
%     for i = 1 : N
% 
%         % Derive the local deterministic cost-of-carry from the integral
%         c_local = (I_ti(i+1) - I_ti(i)) / h;
% 
%         for j = 1 : i      % Index over FX nodes
%             for k = 1 : i  % Index over Rate nodes
% 
%                 % Extract the stochastic rate for the current node
%                 x_t = x_tree(i, k);
% 
%                 % Compute the true FX risk-neutral drift rate
%                 drift_rate = x_t + c_local + shift_amount;
% 
%                 % Calculate the exact risk-neutral expected future value (continuous)
%                 drift = drift_rate * S_tree_rob(i, j);
%                 threshold = S_tree_rob(i, j) + drift * h;
% 
%                 % Continuous vs Discrete Mapping 
%                 % Search the t+1 grid to find the immediate lower floor (jd) 
%                 % and immediate upper ceiling (ju) nodes.
%                 jd = find(S_tree_rob(i+1, 1:i+1) <= threshold, 1, 'last');
%                 ju = find(S_tree_rob(i+1, 1:i+1) >= threshold, 1, 'first');
% 
%                 % Boundary Clamping (Tail Truncation) 
%                 % If an extreme interest rate shock pushes the expected FX value 
%                 % completely outside the physical boundaries, clamp to outermost nodes.
%                 if isempty(jd), jd = 1; end
%                 if isempty(ju), ju = i+1; end
% 
%                 jd_S_rob(i, j, k) = jd;
%                 ju_S_rob(i, j, k) = ju;
% 
%                 % The Collapse Failsafe (ju = jd) 
%                 if ju > jd
%                     % Standard linear interpolation to center the expected value
%                     num = threshold - S_tree_rob(i+1, jd);
%                     den = S_tree_rob(i+1, ju) - S_tree_rob(i+1, jd);
%                     p_hat_S_rob(i, j, k) = max(0, min(1, num/den));
%                 else
%                     % If the ceiling and floor collapse into the exact same node
%                     p_hat_S_rob(i, j, k) = 0.5;
%                 end
% 
%             end
%         end
%     end
% end