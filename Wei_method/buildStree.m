% function [S_tree, jd_S, ju_S, p_hat_S] = buildStree(N, h, kappa, sigma_S, ...
%                    sigma_r, S0, Y0, rho, I, zeta_tree, yearFrac, shift_amount)
% % BUILDSTREE Constructs the bivariate binomial lattice for the FX underlying
% % using the Wei orthogonalization method.
% % It builds an auxiliary grid (Y) driven by the pure stochastic rate (zeta),
% % and adds the deterministic market forward curve (I) at the end.
% %
% % INPUTS:
% %   N            - Number of time steps for the binomial tree discretization.
% %   h            - Length of each time interval (dt = T/N).
% %   kappa        - Mean reversion speed parameter of the Hull-White short rate process.
% %   sigma_S      - Constant volatility of the Foreign Exchange (FX) rate process.
% %   sigma_r      - Constant volatility of the domestic stochastic short rate process.
% %   S0           - Initial spot FX rate (e.g., USD per 1 EUR).
% %   Y0           - Initial value of the auxiliary orthogonalized process Y.
% %   rho           - Correlation coefficient between the FX and short rate Brownian motions.
% %   I            - Vector containing the pre-calculated deterministic forward curve integration components.
% %   zeta_tree    - 2D matrix representing the discretized lattice of the transformed short rate process.
% %   yearFrac     - Time structure vector associated with the market curve data points.
% %   shift_amount - Optional scalar for yield curve parallel shifts.
% %
% % OUTPUTS:
% %   S_tree       - 3D tensor [time, FX-node, short-rate-node] of reconstructed FX spot prices.
% %   jd_S         - 3D tensor of lattice indices mapping the downward transitions on the Y-grid.
% %   ju_S         - 3D tensor of lattice indices mapping the upward transitions on the Y-grid.
% %   p_hat_S      - 3D tensor storing the risk-neutral transition probabilities for the FX component.
% 
%     if nargin < 12
%         shift_amount = 0;
%     end
% 
%     % --- 1. Auxiliary Y-Tree Initialization ---
%     % Rigid, symmetric grid based only on volatility (unit diffusion)
%     Y_tree = NaN(N+1, N+1);
%     for i = 1 : N + 1
%         for j = 1 : i
%             Y_tree(i, j) = Y0 + (2*(j - 1) - (i - 1))*sqrt(h);  
%         end
%     end
% 
%     jd_S = NaN(N+1, N+1, N+1);
%     ju_S = NaN(N+1, N+1, N+1);
%     p_hat_S = NaN(N+1, N+1, N+1);
% 
%     % --- 2. Lattice Mapping & Probabilities ---
%     for i = 1 : N
%         for j = 1 : i
%             for k = 1 : i
% 
%                 drift = (1 / sqrt(1-rho^2)) * (sigma_r/sigma_S - rho*kappa) * zeta_tree(i, k);
% 
%                 % Node mapping with boundary failsafes
%                 jd_S(i, j, k) = j + floor((1 + sqrt(h) * drift)/2);
%                 jd_S(i, j, k) = max(1, min(i, jd_S(i, j, k))); 
%                 ju_S(i, j, k) = jd_S(i, j, k) + 1;
% 
%                 % Marginal probabilities clamped to [0, 1]
%                 num = drift*h + Y_tree(i, j) - Y_tree(i+1, jd_S(i, j, k));
%                 den = 2*sqrt(h);
%                 p_hat_S(i, j, k) = max(0, min(1, num/den));
%             end
%         end
%     end
% 
%     S_tree = NaN(N+1, N+1, N+1);
%     t_nodes = (0:N) * h;
% 
%     % Interpolate the market integral I(t) to match lattice time-steps
%     I_ti = interp1(yearFrac, I, t_nodes, 'linear', 'extrap');
%     I_ti = I_ti + (shift_amount * t_nodes);
% 
%     % --- 4. Final FX Spot Tree Reconstruction ---
%     % Recombine Y-grid, stochastic rate (zeta), and market curve (I_ti)
%     for i = 1 : N + 1 
%         for j = 1 : i
%             for k = 1 : i
%                 S_tree(i, j, k) = S0 * exp(sigma_S * (sqrt(1-rho^2) * Y_tree(i, j) ...
%                     + rho * zeta_tree(i, k)) + I_ti(i) - (sigma_S^2 / 2) * t_nodes(i));
%             end
%         end
%     end
% end

%% Vectorialized function

function [S_tree, jd_S, ju_S, p_hat_S] = buildStree(N, h, kappa, sigma_S, ...
                   sigma_r, S0, Y0, rho, I, zeta_tree, yearFrac, shift_amount)
% BUILDSTREE Constructs the bivariate binomial lattice for the FX underlying
% using the Wei orthogonalization method.
% It builds an auxiliary grid (Y) driven by the pure stochastic rate (zeta),
% and adds the deterministic market forward curve (I) at the end.

    if nargin < 12
        shift_amount = 0;
    end
    
    % --- 1. Auxiliary Y-Tree Initialization VECTORIZED ---
    % Rigid, symmetric grid based only on volatility (unit diffusion)
    % We replace the i and j loops with 2D grid generation
    i_grid = (1:N+1)'; % Column vector representing time indices
    j_grid = 1:N+1;    % Row vector representing space indices
    
    % Compute the whole N+1 x N+1 matrix simultaneously using implicit expansion
    Y_tree_full = Y0 + (2*(j_grid - 1) - (i_grid - 1)) * sqrt(h);  
    
    % Mask out the upper triangle to keep only the valid binomial nodes (where j <= i)
    Y_tree = NaN(N+1, N+1);
    valid_nodes_mask = (j_grid <= i_grid);
    Y_tree(valid_nodes_mask) = Y_tree_full(valid_nodes_mask);
    
    % Preallocate 3D output tensors
    jd_S = NaN(N+1, N+1, N+1);
    ju_S = NaN(N+1, N+1, N+1);
    p_hat_S = NaN(N+1, N+1, N+1);
    
    % Pre-calculate the constant scalar portion of the drift for efficiency
    drift_scalar = (1 / sqrt(1-rho^2)) * (sigma_r/sigma_S - rho*kappa);
    
    % --- 2. Lattice Mapping & Probabilities VECTORIZED ---
    % We keep only the time loop 'i', calculating entire (i x i) spatial planes at once
    for i = 1 : N
        
        % Create orthogonal 1D vectors for indices j and k
        j_vec = (1:i)';                 % Column vector (i x 1) -> represents FX dimension
        zeta_vec = zeta_tree(i, 1:i);   % Row vector (1 x i)    -> represents Rate dimension
        
        % Vectorized drift: depends only on 'k', so it remains a row vector (1 x i)
        drift_vec = drift_scalar * zeta_vec;
        
        % Node mapping via Implicit Expansion: 
        % (Column i x 1) + (Row 1 x i) automatically generates a full (i x i) matrix
        jd_mat = j_vec + floor((1 + sqrt(h) * drift_vec)/2);
        
        % Apply boundary failsafes simultaneously to the entire matrix
        jd_mat = max(1, min(i, jd_mat)); 
        ju_mat = jd_mat + 1;
        
        % --- Probability Computation ---
        % Extract current Y nodes (depends on j)
        Y_cur = Y_tree(i, 1:i)'; % Column vector (i x 1)
        
        % Extract future Y nodes:
        % We take the 1D row from time i+1, and index it using our 2D jd_mat.
        % MATLAB handles this perfectly, returning an (i x i) matrix of future Y values!
        Y_next_row = Y_tree(i+1, 1:i+1); 
        Y_next_mat = Y_next_row(jd_mat); % Result is (i x i)
        
        % Compute numerator: (Row 1xi) + (Col ix1) - (Matrix ixi) = Matrix ixi
        num = drift_vec*h + Y_cur - Y_next_mat;
        den = 2*sqrt(h);
        
        % Marginal probabilities clamped to [0, 1]
        p_hat_mat = max(0, min(1, num/den));
        
        % Inject the computed 2D planes into the 3D tensors using reshape
        jd_S(i, 1:i, 1:i) = reshape(jd_mat, [1, i, i]);
        ju_S(i, 1:i, 1:i) = reshape(ju_mat, [1, i, i]);
        p_hat_S(i, 1:i, 1:i) = reshape(p_hat_mat, [1, i, i]);
    end
    
    % --- 3. Market Integral Interpolation ---
    S_tree = NaN(N+1, N+1, N+1);
    t_nodes = (0:N) * h;
    
    % Interpolate the market integral I(t) to match lattice time-steps
    I_ti = interp1(yearFrac, I, t_nodes, 'linear', 'extrap');
    I_ti = I_ti + (shift_amount * t_nodes);
    
    % --- 4. Final FX Spot Tree Reconstruction VECTORIZED ---
    % Recombine Y-grid, stochastic rate (zeta), and market curve (I_ti)
    % We process the entire (i x i) spatial slice for each time step i
    for i = 1 : N + 1 
        
        % Extract 1D vectors
        Y_cur = Y_tree(i, 1:i)';      % Column vector (i x 1)
        zeta_cur = zeta_tree(i, 1:i); % Row vector (1 x i)
        
        % The exponential argument instantly becomes an (i x i) matrix
        exponent = sigma_S * (sqrt(1-rho^2) * Y_cur + rho * zeta_cur) ...
                   + I_ti(i) - (sigma_S^2 / 2) * t_nodes(i);
               
        % Calculate the final spot prices for the current slice
        S_slice = S0 * exp(exponent);
        
        % Insert the matrix into the 3D tensor
        S_tree(i, 1:i, 1:i) = reshape(S_slice, [1, i, i]);
    end
end