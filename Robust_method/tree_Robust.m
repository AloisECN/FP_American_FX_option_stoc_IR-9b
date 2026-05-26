function [S_tree, q_ju_ku, q_ju_kd, q_jd_ku, q_jd_kd, ju_S, jd_S, ku_r, kd_r, df_rob] = ...
    tree_Robust(N, h, kappa, zeta0, sigma_S, sigma_r, S0, rho, I, yearFrac)

%% Part 1: Interest Rate Tree
r_tree = zeros(N+1, N+1);
for i = 0:N
    for k = 0:i
        %original zeta tree independant from S
        val = zeta0 + (2*k - i) * sqrt(h);
        if val >= 0
            r_tree(i+1, k+1) = (val * sigma_r)^2 / 4;
        else
            r_tree(i+1, k+1) = 0;
        end
    end
end

% Compute discount factors
df_rob = exp(-r_tree * h);   % (N+1)x(N+1), same shape as r_tree

kd_r = NaN(N+1, N+1);
ku_r = NaN(N+1, N+1);
p_r  = NaN(N+1, N+1);

% Building probabilities equation 13 , 14 , 15 from Apploloni 

for i = 0:N-1
    for k = 0:i
        drift     = -kappa * r_tree(i+1, k+1); % adapted to our HW process
        threshold = r_tree(i+1, k+1) + drift * h;

        %Eq 13
        cands_d = 0:k;
        mask_d  = r_tree(i+2, cands_d + 1) <= threshold;
        if any(mask_d)
            kd_r(i+1, k+1) = max(cands_d(mask_d));
        else
            kd_r(i+1, k+1) = 0;
        end

        %Eq 14
        cands_u = (k+1):(i+1);
        mask_u  = r_tree(i+2, cands_u + 1) >= threshold;
        if any(mask_u)
            ku_r(i+1, k+1) = min(cands_u(mask_u));
        else
            ku_r(i+1, k+1) = i+1;
        end

        %Eq 15
        num = drift * h + r_tree(i+1, k+1) - r_tree(i+2, kd_r(i+1,k+1) + 1);
        den = r_tree(i+2, ku_r(i+1,k+1) + 1) - r_tree(i+2, kd_r(i+1,k+1) + 1);
        p_r(i+1, k+1) = max(0, min(1, num/den));

    end
end

%% Part 2: FX Tree 
% Eq 16 and 17 
U0     = log(S0) / sigma_S;
U_tree = zeros(N+1, N+1);
for i = 0:N
    for j = 0:i
        U_tree(i+1, j+1) = U0 + (2*j - i) * sqrt(h);
    end
end
S_tree = exp(sigma_S * U_tree);

%Computing probabilities

jd_S = NaN(N+1, N+1, N+1);
ju_S = NaN(N+1, N+1, N+1);
p_S  = NaN(N+1, N+1, N+1);

for i = 0:N-1
    for j = 0:i
        for k = 0:i
            S_cur     = S_tree(i+1, j+1);
            r_cur     = r_tree(i+1, k+1);
            drift     = S_cur * r_cur; % drift change => normal as r_t changes too
            threshold = S_cur + drift * h;

            %Eq 18
            cands_d = 0:j;
            mask_d  = S_tree(i+2, cands_d + 1) <= threshold;
            if any(mask_d)
                jd_S(i+1, j+1, k+1) = max(cands_d(mask_d));
            else
                jd_S(i+1, j+1, k+1) = 0;
            end

            %Eq 19
            cands_u = (j+1):(i+1);
            mask_u  = S_tree(i+2, cands_u + 1) >= threshold;
            if any(mask_u)
                ju_S(i+1, j+1, k+1) = min(cands_u(mask_u));
            else
                ju_S(i+1, j+1, k+1) = i+1;
            end

            %Eq 20
            num = drift * h + S_cur - S_tree(i+2, jd_S(i+1,j+1,k+1) + 1);
            den = S_tree(i+2, ju_S(i+1,j+1,k+1) + 1) - S_tree(i+2, jd_S(i+1,j+1,k+1) + 1);
            p_S(i+1, j+1, k+1) = max(0, min(1, num/den));
        
        end
    end
end

%% Part 3: Joint probabilities

q_ju_ku = NaN(N+1, N+1, N+1);
q_ju_kd = NaN(N+1, N+1, N+1);
q_jd_ku = NaN(N+1, N+1, N+1);
q_jd_kd = NaN(N+1, N+1, N+1);

for i = 0:N-1
    for j = 0:i
        for k = 0:i
            %Declare all varibales for clarity 
            %Not memory fficient but who cares 

            ju_t  = ju_S(i+1, j+1, k+1);
            jd_t  = jd_S(i+1, j+1, k+1);
            ku_t  = ku_r(i+1, k+1);
            kd_t  = kd_r(i+1, k+1);

            S_cur = S_tree(i+1, j+1);
            r_cur = r_tree(i+1, k+1);

            S_ju = S_tree(i+2, ju_t + 1);
            S_jd = S_tree(i+2, jd_t + 1);
            r_ku = r_tree(i+2, ku_t + 1);
            r_kd = r_tree(i+2, kd_t + 1);

            %Eq 22
            m_ju_ku = (S_ju - S_cur) * (r_ku - r_cur);
            m_ju_kd = (S_ju - S_cur) * (r_kd - r_cur);
            m_jd_ku = (S_jd - S_cur) * (r_ku - r_cur);
            m_jd_kd = (S_jd - S_cur) * (r_kd - r_cur);
            
            %Eq 21
            %RHS of (21)
            p_hat_i_j_k = p_S(i+1, j+1, k+1); % first
            p_i_k     = p_r(i+1, k+1); %second
            C = rho * sigma_r *  sqrt(r_cur) * sigma_S * S_cur * h ; %fourth
            
            %Solving the system explicitly yields
            
            A = [ 1,1,0,0; ...
                1,0,1,0; ...
                1,1,1,1; ...
                m_ju_ku, m_ju_kd, m_jd_ku, m_jd_kd ];

            % det(A) = -m_{ju\_ku} + m_{ju\_kd} + m_{jd\_ku} - m_{jd\_kd)
            
            B = [ p_hat_i_j_k ; p_i_k ; 1 ; C ];

            q_vec = lsqnonneg(A, B);

            q_ju_ku(i+1,j+1,k+1) = q_vec(1);
            q_ju_kd(i+1,j+1,k+1) = q_vec(2);
            q_jd_ku(i+1,j+1,k+1) = q_vec(3);
            q_jd_kd(i+1,j+1,k+1) = q_vec(4);
        end
    end
end

end 