function [S_tree, q_ju_ku, q_ju_kd, q_jd_ku, q_jd_kd, ju_S, jd_S, ku_r, kd_r, df_rob] = ...
    tree_Robust(N, h, kappa, zeta0, sigma_S, sigma_r, S0, rho, I, yearFrac)

%% Part 1: Interest Rate Tree
r_tree = zeros(N+1, N+1);
for i = 0:N
    for k = 0:i
        val = zeta0 + (2*k - i) * sqrt(h);
        if val >= 0
            r_tree(i+1, k+1) = (val * sigma_r)^2 / 4;
        else
            r_tree(i+1, k+1) = 0;
        end
    end
end

% --- SELF-CONSISTENT DISCOUNT FACTORS for the Robust r_tree ---
% df_rob(i,k) = exp(-r_tree(i,k)*h), used in backward induction at node (i,k)
df_rob = exp(-r_tree * h);   % (N+1)x(N+1), same shape as r_tree

kd_r = NaN(N+1, N+1);
ku_r = NaN(N+1, N+1);
p_r  = NaN(N+1, N+1);

for i = 0:N-1
    for k = 0:i
        drift     = -kappa * r_tree(i+1, k+1);
        threshold = r_tree(i+1, k+1) + drift * h;

        cands_d = 0:k;
        mask_d  = r_tree(i+2, cands_d + 1) <= threshold;
        if any(mask_d)
            kd_r(i+1, k+1) = max(cands_d(mask_d));
        else
            kd_r(i+1, k+1) = 0;
        end

        cands_u = (k+1):(i+1);
        mask_u  = r_tree(i+2, cands_u + 1) >= threshold;
        if any(mask_u)
            ku_r(i+1, k+1) = min(cands_u(mask_u));
        else
            ku_r(i+1, k+1) = i+1;
        end

        num = drift * h + r_tree(i+1, k+1) - r_tree(i+2, kd_r(i+1,k+1) + 1);
        den = r_tree(i+2, ku_r(i+1,k+1) + 1) - r_tree(i+2, kd_r(i+1,k+1) + 1);
        if den ~= 0
            p_r(i+1, k+1) = max(0, min(1, num/den));
        else
            p_r(i+1, k+1) = 0.5;
        end
    end
end

%% Part 2: FX Tree
U0     = log(S0) / sigma_S;
U_tree = zeros(N+1, N+1);
for i = 0:N
    for j = 0:i
        U_tree(i+1, j+1) = U0 + (2*j - i) * sqrt(h);
    end
end
S_tree = exp(sigma_S * U_tree);

jd_S = NaN(N+1, N+1, N+1);
ju_S = NaN(N+1, N+1, N+1);
p_S  = NaN(N+1, N+1, N+1);

for i = 0:N-1
    for j = 0:i
        for k = 0:i
            S_cur     = S_tree(i+1, j+1);
            r_cur     = r_tree(i+1, k+1);
            drift     = S_cur * r_cur;
            threshold = S_cur + drift * h;

            cands_d = 0:j;
            mask_d  = S_tree(i+2, cands_d + 1) <= threshold;
            if any(mask_d)
                jd_S(i+1, j+1, k+1) = max(cands_d(mask_d));
            else
                jd_S(i+1, j+1, k+1) = 0;
            end

            cands_u = (j+1):(i+1);
            mask_u  = S_tree(i+2, cands_u + 1) >= threshold;
            if any(mask_u)
                ju_S(i+1, j+1, k+1) = min(cands_u(mask_u));
            else
                ju_S(i+1, j+1, k+1) = i+1;
            end

            num = drift * h + S_cur - S_tree(i+2, jd_S(i+1,j+1,k+1) + 1);
            den = S_tree(i+2, ju_S(i+1,j+1,k+1) + 1) - S_tree(i+2, jd_S(i+1,j+1,k+1) + 1);
            if den ~= 0
                p_S(i+1, j+1, k+1) = max(0, min(1, num/den));
            else
                p_S(i+1, j+1, k+1) = 0.5;
            end
        end
    end
end

%% Part 3: Joint probabilities with ENFORCED normalization
q_ju_ku = NaN(N+1, N+1, N+1);
q_ju_kd = NaN(N+1, N+1, N+1);
q_jd_ku = NaN(N+1, N+1, N+1);
q_jd_kd = NaN(N+1, N+1, N+1);

for i = 0:N-1
    for j = 0:i
        for k = 0:i
            ju_t  = ju_S(i+1, j+1, k+1);
            jd_t  = jd_S(i+1, j+1, k+1);
            ku_t  = ku_r(i+1, k+1);
            kd_t  = kd_r(i+1, k+1);

            p_hat = p_S(i+1, j+1, k+1);
            p     = p_r(i+1, k+1);

            S_cur = S_tree(i+1, j+1);
            r_cur = r_tree(i+1, k+1);

            S_ju = S_tree(i+2, ju_t + 1);
            S_jd = S_tree(i+2, jd_t + 1);
            r_ku = r_tree(i+2, ku_t + 1);
            r_kd = r_tree(i+2, kd_t + 1);

            m_ju_ku = (S_ju - S_cur) * (r_ku - r_cur);
            m_ju_kd = (S_ju - S_cur) * (r_kd - r_cur);
            m_jd_ku = (S_jd - S_cur) * (r_ku - r_cur);
            m_jd_kd = (S_jd - S_cur) * (r_kd - r_cur);

            C     = rho * sigma_r * sigma_S * S_cur * h;
            denom = m_ju_ku - m_ju_kd - m_jd_ku + m_jd_kd;
            numer = C - m_ju_kd*p_hat - m_jd_ku*p - m_jd_kd*(1 - p_hat - p);

            if abs(denom) > 1e-14
                a = numer / denom;
            else
                a = p_hat * p;  % independence fallback
            end

            % Recover the four joint probabilities
            a = p_hat * p;     % q(ju, ku)  — fallback to independence when covariance solve is ill-conditioned
            b = p_hat - a;     % q(ju, kd)
            c = p     - a;     % q(jd, ku)
            d = 1 - p_hat - p + a; % q(jd, kd)

            % --- CRITICAL: floor at 0, then RENORMALIZE to enforce sum = 1 ---
            q_raw = max(0, [a, b, c, d]);
            s = sum(q_raw);
            if s > 1e-14
                q_raw = q_raw / s;
            else
                q_raw = [p_hat*p, p_hat*(1-p), (1-p_hat)*p, (1-p_hat)*(1-p)];
            end

            q_ju_ku(i+1,j+1,k+1) = q_raw(1);
            q_ju_kd(i+1,j+1,k+1) = q_raw(2);
            q_jd_ku(i+1,j+1,k+1) = q_raw(3);
            q_jd_kd(i+1,j+1,k+1) = q_raw(4);
        end
    end
end

end