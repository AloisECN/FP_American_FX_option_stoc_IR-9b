%% Vectorialized function

function [x_tree, zeta_tree, kd_r, ku_r, p_r] = buildXtree(N, h, kappa, sigma_r, zeta0)
% BUILDXTREE Constructs the discrete binomial lattice for the stochastic 
% component of the short-rate process (e.g., Hull-White model).
%
% This function first builds a normalized, dimensionless tree (zeta) with unit 
% variance. It then calculates the mean-reverting drift, the branching indices 
% (jumps), and the transition probabilities. Finally, it scales the tree by 
% the volatility to obtain the actual stochastic rate process (x_tree).
%
% INPUTS:
%   N       - Number of time steps in the tree.
%   h       - Time step size (dt = T/N).
%   kappa   - Mean reversion speed of the short-rate process.
%   sigma_r - Constant volatility of the short-rate process.
%   zeta0   - Initial value of the normalized process (usually 0).
%
% OUTPUTS:
%   x_tree    - 2D matrix representing the stochastic rate process (scaled by sigma_r).
%   zeta_tree - 2D matrix representing the normalized, dimensionless grid.
%   kd_r      - 2D matrix of indices mapping the downward transitions on the tree.
%   ku_r      - 2D matrix of indices mapping the upward transitions on the tree.
%   p_r       - 2D matrix storing the risk-neutral upward transition probabilities.

    %  1. Normalized Grid (zeta_tree) Initialization VECTORIZED 
    % We replace the nested loops with MATLAB's implicit expansion.
    % We create orthogonal vectors for time (i) and space (k).
    i_grid = (1:N+1)';  % Column vector for time dimension
    k_grid = 1:N+1;     % Row vector for space dimension
    
    % Compute the full matrix instantly
    zeta_full = zeta0 + (2*k_grid - i_grid - 1) * sqrt(h);
    
    % Apply a logical mask to keep only the lower triangular part (where k <= i)
    % which represents the valid nodes of our binomial tree.
    zeta_tree = zeros(N+1, N+1);
    valid_nodes_mask = (k_grid <= i_grid);
    zeta_tree(valid_nodes_mask) = zeta_full(valid_nodes_mask);

    %  2. Preallocate Arrays for Jumps and Probabilities 
    kd_r = NaN(N+1, N+1);
    ku_r = NaN(N+1, N+1);
    p_r  = NaN(N+1, N+1);

    %  3. Lattice Mapping & Probabilities VECTORIZED 
    % We eliminate the inner 'k' loop. We process an entire time step 'i' 
    % (a 1D array of nodes) in a single vectorized operation.
    for i = 1 : N
        
        % Extract current nodes and generate space indices vector
        k_vec = 1:i;                    % Row vector of indices 1 to i
        zeta_cur = zeta_tree(i, 1:i);   % Row vector of zeta values at time i
        
        % Compute mean-reverting drift for all nodes at time i simultaneously
        drift = -kappa * zeta_cur;
        
        % Calculate downward jump indices (kd) for the whole vector
        kd_vec = k_vec + floor((drift * sqrt(h) + 1) / 2);

        % Apply boundary conditions (capping) to ensure indices remain within the tree.
        % 1 is the lowest index, i is the maximum possible downward jump boundary.
        kd_vec = max(1, min(i, kd_vec));
        
        % Upward jump is deterministically the node immediately above the downward jump
        ku_vec = kd_vec + 1;

        %  Probability Computation 
        % Extract the whole row of zeta values at the NEXT time step (i+1)
        zeta_next = zeta_tree(i+1, 1:i+1);
        
        % Vectorized indexing: fetch the future zeta values corresponding to our downward jumps
        zeta_kd = zeta_next(kd_vec);

        % Compute probabilities using the mathematical formula
        num = drift * h + zeta_cur - zeta_kd;
        den = 2 * sqrt(h);
        
        % Clamp probabilities to the valid [0, 1] range to avoid negative probabilities
        p_r_vec = max(0, min(1, num ./ den));
        
        % Store the vectorized results back into the matrices
        kd_r(i, 1:i) = kd_vec;
        ku_r(i, 1:i) = ku_vec;
        p_r(i, 1:i)  = p_r_vec;
    end

    %  4. Final Stochastic Process (x_tree) 
    % Scale the normalized tree by the volatility to obtain the actual x_t process.
    % This is a simple scalar-matrix multiplication.
    x_tree = sigma_r * zeta_tree;

end

% Legacy 

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
% 
% 
