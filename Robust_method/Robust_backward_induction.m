function V0 = Robust_backward_Induction(S_tree, N, K, ...
                 q_ju_ku, q_ju_kd, q_jd_ku, q_jd_kd, ...
                 ju_S, jd_S, ku_r, kd_r, df_rob)

    V = NaN(N+1, N+1, N+1);
    for j = 1:N+1
        for k = 1:N+1
            V(N+1, j, k) = max(K - S_tree(N+1, j), 0);
        end
    end

    %backward induction loop
    for i = N:-1:1
        for j = 1:i
            for k = 1:i
                p_uu = q_ju_ku(i, j, k);
                p_ud = q_ju_kd(i, j, k);
                p_du = q_jd_ku(i, j, k);
                p_dd = q_jd_kd(i, j, k);

                ju = ju_S(i, j, k) + 1;
                jd = jd_S(i, j, k) + 1;
                ku = ku_r(i, k)    + 1;
                kd = kd_r(i, k)    + 1;

                % df_rob(i,k): discount from node (i,k) to next step
                df = df_rob(i, k);
                continuation = df * ( p_uu * V(i+1, ju, ku) ...
                                    + p_ud * V(i+1, ju, kd) ...
                                    + p_du * V(i+1, jd, ku) ...
                                    + p_dd * V(i+1, jd, kd) );

                intrinsic  = max(K - S_tree(i, j), 0);
                V(i, j, k) = max(intrinsic, continuation); % as it's an american option
            end
        end
    end

    V0 = V(1, 1, 1);
end