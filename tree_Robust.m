function price  = tree_Robust(N, h, kappa,zeta0, sigma_S, sigma_r, S0, X0, ro, I, yearFrac)

% Implement tree with the robust method
% inputs
%
%
% Ouputs
%   Price of the Put option

% First part modified r_tree tree

r_tree = zeros(N+1, N+1);
for i = 0:N
    for k = 0:i
        if (zeta0 + (2*k - i)) >= 0 
            r_tree(i+1, k+1) = (zeta0 + (2*k - i) * sqrt(h))**2 * sigma_r**2 / 4 ;
        else
            r_tree(i+1, k+1) = 0 ;
    end
end

kd_r = NaN(N+1, N+1);
ku_r = NaN(N+1, N+1);
p_r  = NaN(N+1, N+1);

for i = 0:N-1
    for k = 0:i

        drift = -kappa * r_tree(i+1, k+1);

        threshold = r_tree(i+1, k+1) + drift * h;

        k_star_candidates = 0 : k;           
        mask_d = r_tree(i+1, k_star_candidates + 1) <= threshold; 
        if any(mask_d)
            kd(i+1, k+1) = max(k_star_candidates(mask_d));
        else
            kd(i+1, k+1) = NaN;  % no valid k* found
        end

        k_star_candidates = (k + 1) : (i + 1);
        mask_u = r_tree(i+1, k_star_candidates + 1) >= threshold; 

        if any(mask_u)
            ku(i+1, k+1) = min(k_star_candidates(mask_u));
        else
            ku(i+1, k+1) = NaN;  % no valid k* found
        end
         
        
        num = drift *  h + r_tree (i+1, k + 1) - r_tree(i +1 , kd_r(i+1, k+1))
        den =  r_tree(i +1 , ku_r(i+1, k+1)) - r_tree(i +1 , kd_r(i+1, k+1))
        p_r(i+1, k+1) = max(0, min(1, num/den));
    end

end

%Second Part the 2D tree 
%simple U-tree

U0 = 1/sigma_S * log(S0)
U_tree = zeros(N+1, N+1);
for i = 0:N
    for j = 0:i
        U_tree[i+1, j+1] = U0 + (2*j - i) * sqrt(h)
    end
end

S_tree = exp(sigma_S * U_tree)

% Compute as independant 

jd = NaN(N+1,N+1, N+1);
ju = NaN(N+1,N+1, N+1);
p_S  = NaN(N+1,N+1, N+1);

for i = 0:N-1
    for j = 0:i
        for k = 0:i
            drift = 


            num = 
            den = 
            p_r(i+1, j+1, k+1) = max(0, min(1, num/den));
        end
    end 
end

% Introduce the covariance struct 
qd = NaN(N+1,N+1, N+1);
qu = NaN(N+1,N+1, N+1);
p_S  = NaN(N+1,N+1, N+1);

end 