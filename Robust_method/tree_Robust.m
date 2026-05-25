function []  = tree_Robust(N, h, kappa,zeta0, sigma_S, sigma_r, S0, X0, ro, I, yearFrac)

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

            drift = S_tree(i,j) * r_tree(i, k);
            threshold = S_tree(i,j) + drift * h;
            j_star_candidates = 0 : j;           
            mask_d = S_tree(i+1, j_star_candidates + 1) <= threshold; 
            if any(mask_d)
                jd(i+1,j+1, k+1) = max(j_star_candidates(mask_d));
            else
                jd(i+1,j+1, k+1) = NaN;  % no valid k* found
            end

            j_star_candidates = (j + 1) : (i + 1);
            mask_u = S_tree(i+1, j_star_candidates + 1) >= threshold; 

            if any(mask_u)
                ju(i+1,j+1, k+1) = min(j_star_candidates(mask_u));
            else
                ju(i+1,j+1, k+1) = NaN;  % no valid k* found
            end

            num = drift * h + S_tree(i,j) - S_tree(i+1 , jd(i+1, j+1 , k+1));

            den = S_tree(i+1 , ju(i+1, j+1 , k+1)) - S_tree(i+1 , jd(i+1, j+1 , k+1));
             
            p_S(i+1, j+1, k+1) = max(0, min(1, num/den));
        end
    end 
end

% Introduce the covariance struct 

% Transition probabilities for the bivariate tree (covariance matching)
q_ju_ku = NaN(N+1, N+1, N+1);
q_ju_kd = NaN(N+1, N+1, N+1);
q_jd_ku = NaN(N+1, N+1, N+1);
q_jd_kd = NaN(N+1, N+1, N+1);

for i = 0:N-1
    for j = 0:i
        for k = 0:i

            % Recover proba from prev trees
            ju_temp = ju(i+1, j+1, k+1);
            jd_temp = jd(i+1, j+1, k+1);
            ku_temp = ku_r(i+1, k+1);
            kd_temp = kd_r(i+1, k+1);

            p_hat = p_S(i+1, j+1, k+1);  
            p     = p_r(i+1, k+1);       

            % Current node values
            S_cur = S_tree(i+1, j+1);
            r_cur = r_tree(i+1, k+1);

            % some allocation for clarity
            S_ju = S_tree(i+2, ju_temp + 1);
            S_jd = S_tree(i+2, jd_temp + 1);
            r_ku = r_tree(i+2, ku_temp + 1);
            r_kd = r_tree(i+2, kd_temp + 1);

            % The four "m" terms from eq. (22)
            m_ju_ku = (S_ju - S_cur) * (r_ku - r_cur);
            m_ju_kd = (S_ju - S_cur) * (r_kd - r_cur);
            m_jd_ku = (S_jd - S_cur) * (r_ku - r_cur);
            m_jd_kd = (S_jd - S_cur) * (r_kd - r_cur);

            % RHS of eq. 4 in system (21)
            C = ro * sigma_r * sigma_S * S_cur * h;

            % Solve for q_ju_ku from the 4th equation
            
            denom = m_ju_ku - m_ju_kd - m_jd_ku + m_jd_kd;
            numer = C - m_ju_kd * p_hat - m_jd_ku * p - m_jd_kd * (1 - p_hat - p);

            a = numer / denom;

            % Recover the other three probabilities
            b = p_hat - a;   % q(ju, kd)
            c = p     - a;   % q(jd, ku)
            d = 1 - p_hat - p + a;  % q(jd, kd)

            q_ju_ku(i+1, j+1, k+1) = max(0, min(1, a));
            q_ju_kd(i+1, j+1, k+1) = max(0, min(1, b));
            q_jd_ku(i+1, j+1, k+1) = max(0, min(1, c));
            q_jd_kd(i+1, j+1, k+1) = max(0, min(1, d));

        end
    end
end

end % of the fuction 