function V0 = Robust_backward_Induction(S_tree, N, K, ...
                 q_ju_ku, q_ju_kd, q_jd_ku, q_jd_kd, ...
                 ju_S, jd_S, ku_r, kd_r, df_HW)

    % Terminal payoff at maturity (i = N, so layer N+1 in 1-based)
    V = NaN(N+1, N+1, N+1);
    for j = 1 : N+1
        for k = 1 : N+1
            V(N+1, j, k) = max(K - S_tree(N+1, j), 0);
            % Note: S_tree is (i,j) not (i,j,k) — S doesn't depend on k
        end
    end

    % Backward induction
    for i = N : -1 : 1
        for j = 1 : i
            for k = 1 : i

                % Joint probabilities from system (21)
                p_uu = q_ju_ku(i, j, k);
                p_ud = q_ju_kd(i, j, k);
                p_du = q_jd_ku(i, j, k);
                p_dd = q_jd_kd(i, j, k);

                % Next-step indices (already 1-based if stored that way,
                % or add +1 here if stored as 0-based)
                ju = ju_S(i, j, k) + 1;
                jd = jd_S(i, j, k) + 1;
                ku = ku_r(i, k)    + 1;
                kd = kd_r(i, k)    + 1;

                % Discounted continuation value
                df = df_HW(i, k);
                continuation = df * ( p_uu * V(i+1, ju, ku) ...
                                    + p_ud * V(i+1, ju, kd) ...
                                    + p_du * V(i+1, jd, ku) ...
                                    + p_dd * V(i+1, jd, kd) );

                % American option: take the better of exercise or hold
                intrinsic = max(K - S_tree(i, j), 0);
                V(i, j, k) = max(intrinsic, continuation);

            end
        end
    end

    % Price at the root node
    V0 = V(1, 1, 1);

end