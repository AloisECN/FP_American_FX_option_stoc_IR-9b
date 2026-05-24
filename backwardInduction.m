function V0 = backwardInduction(S_tree, N, K, p_r, p_hat_S, kd_r, jd_S, df_HW)

V = NaN(N+1, N+1, N+1);
for j = 1 : N+1
    for k = 1: N+1
        V(N+1, j, k) = max(K - S_tree(N+1, j, k), 0);
    end
end

for i = N : -1 : 1
    for j = 1 : i
        for k = 1 : i
            ps = p_hat_S(i, j, k);
            pr = p_r(i, k);

            p_UU = ps * pr;
            p_UD = ps * (1 - pr);
            p_DU = (1 - ps) * pr;
            p_DD = (1 - ps) * (1 - pr);

            jd = jd_S(i, j, k);
            ju = 1 + jd;
            kd = kd_r(i, k);
            ku = 1 + kd;

            V_UU = V(i+1, ju, ku) * p_UU;
            V_UD = V(i+1, ju, kd) * p_UD;
            V_DU = V(i+1, jd, ku) * p_DU;
            V_DD = V(i+1, jd, kd) * p_DD;

            df = df_HW(i, k);
            continuationValue = df * (V_UU + V_UD + V_DU + V_DD);

            intrinsicValue = max(K - S_tree(i, j, k), 0);

            V(i, j, k) = max(intrinsicValue, continuationValue);

        end
    end
end


V0 = V(1, 1, 1);
end