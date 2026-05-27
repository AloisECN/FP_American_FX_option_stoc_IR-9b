function [x_tree, zeta_tree, kd_r, ku_r, p_r] = buildXtree(N, h, kappa, sigma_r, zeta0)

    zeta_tree = zeros(N+1, N+1);
    for i = 1 : N+1
        for k = 1 : i
           zeta_tree(i, k) = zeta0 + (2*k - i - 1) * sqrt(h);
        end
    end

    kd_r = NaN(N+1, N+1);
    ku_r = NaN(N+1, N+1);
    p_r  = NaN(N+1, N+1);

    for i = 1 : N
        for k = 1 : i
            drift = -kappa * zeta_tree(i, k);
            kd_r(i, k) = k + floor((drift * sqrt(h) + 1)/2);

            % 1 is the lowest index of the tree
            kd_r(i, k) = max(1, min(i, kd_r(i, k)));
            ku_r(i,k) = kd_r(i,k) + 1;

            num = drift * h + zeta_tree(i, k) - zeta_tree(i+1, kd_r(i, k));
            den = 2 * sqrt(h);
            p_r(i, k) = max(0, min(1, num/den));
        end

    end

    x_tree = sigma_r * zeta_tree;

end



