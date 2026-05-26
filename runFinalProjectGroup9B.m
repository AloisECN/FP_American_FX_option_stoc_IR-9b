%% MAIN SCRIPT: Project - American FX Option with Stochastic Interest Rate
clc; clear all; close all;
warning('off')

%% 0. MARKET DATA INGESTION & INITIAL CALIBRATION
dates = { '12/Apr/2016'; '13/Apr/2016'; '19/Apr/2016'; '12/May/2016'; ...
          '13/Jun/2016'; '12/Jul/2016'; '12/Oct/2016'; '12/Apr/2017'; ...
          '12/Apr/2018'; '12/Apr/2019'; '14/Apr/2020'};
dates_dt = datetime(dates, 'InputFormat', 'dd/MMM/yyyy', 'Locale', 'en_US');
refDate  = dates_dt(1);
discounts = [ 1; 0.99998157; 0.99987514; 0.99942882; 0.99872386; ...
              0.99787021; 0.99480441; 0.98710348; 0.97399953; ...
              0.95687968; 0.93750149];

fwd_mid_curve = [ 1.14055; 1.14077; 1.14152; 1.14257; 1.14370; ...
                  1.14730; 1.15135; 1.15558; 1.17480; 1.19685 ];
maturities = { '12/Apr/2016'; '19/Apr/2016'; '12/May/2016'; '13/Jun/2016'; ...
               '12/Jul/2016'; '12/Oct/2016'; '12/Jan/2017'; '12/Apr/2017'; ...
               '12/Apr/2018'; '12/Apr/2019'};
maturities_dt = datetime(maturities, 'InputFormat', 'dd/MMM/yyyy', 'Locale', 'en_US');

sigma_r  = 0.01;
kappa    = 0.10;
S0       = 1.14055;
sigma_S  = 0.1048;
rho      = 0.3;

maturity_target = maturities_dt(9);
TTM = yearfrac(refDate, maturity_target, 3);

yearFraction = yearfrac(refDate, maturities_dt, 3);
[costOfCarry, I] = calibrateCostOfCarry(yearFraction, fwd_mid_curve, S0, sigma_r, kappa);

r0    = -1/(1/365) * log(0.99998157);
zeta0 = r0 / sigma_r;

%% Shared parameters (used by both methods and by Point iii)
N     = 100;
h     = TTM / N;
t_nodes  = (0:N)*h;
treeDates = refDate + days(round(365 * t_nodes));

K_atm    = S0;
K_vec    = linspace(0.8 * S0, 1.2 * S0, 100);
rho_new  = linspace(-0.99, 0.99, 50);
N_vec    = [10, 25, 50, 75, 100, 150, 200, 300];

% IR tree and discounts — needed by both methods and by Point iii
B0_nodes = getDiscountFactorByZeroRatesLinearInterp(refDate, treeDates, dates_dt, discounts);
[x_tree, zeta_tree, kd_r, ku_r, p_r] = buildXtree(N, h, kappa, sigma_r, zeta0);
df_HW    = fwdDiscounts_OU(x_tree, kappa, sigma_r, t_nodes, B0_nodes);
X0       = -rho * zeta0 / sqrt(1 - rho^2);

%% Method selection
method = input('Which method? 1: Wei  |  2: Robust  --> ');

if method == 1

    %% POINT ii-a: American Put Price vs Strike (Wei)
    [S_tree, jd_S, ju_S, p_hat_S] = buildStree(N, h, kappa, sigma_S, sigma_r, ...
                                                S0, X0, rho, I, zeta_tree, yearFraction);

    prices_vs_K = zeros(size(K_vec));
    for idx = 1:length(K_vec)
        prices_vs_K(idx) = backwardInduction(S_tree, N, K_vec(idx), p_r, p_hat_S, kd_r, jd_S, df_HW);
    end

    figure('Name', 'American Put vs Strike');
    plot(K_vec, prices_vs_K, '-o', 'LineWidth', 1.5); grid on;
    title('American Put Option Price vs Strike (T=2y, \rho=0.3)');
    xlabel('Strike Price (K)'); ylabel('Option Price (USD)');
    xline(S0, '--r', 'ATM Strike');
    legend('American Put Price', 'ATM Level', 'Location', 'best');

    %% POINT ii-b: ATM Price vs Correlation (Wei)
    prices_vs_rho = zeros(size(rho_new));
    for i = 1:length(rho_new)
        current_rho   = rho_new(i);
        X0_current    = -current_rho * zeta0 / sqrt(1 - current_rho^2);
        [S_tree_curr, jd_S_curr, ju_S_curr, p_hat_S_curr] = buildStree(N, h, kappa, sigma_S, sigma_r, ...
                                                                        S0, X0_current, current_rho, I, zeta_tree, yearFraction);
        prices_vs_rho(i) = backwardInduction(S_tree_curr, N, K_atm, p_r, p_hat_S_curr, kd_r, jd_S_curr, df_HW);
    end

    figure('Name', 'American Put vs Correlation');
    plot(rho_new, prices_vs_rho, '-ro', 'LineWidth', 1.5); grid on;
    title('ATM American Put Price vs Correlation \rho (T=2y, K=S_0)');
    xlabel('Correlation (\rho)'); ylabel('Option Price (USD)');

    %% POINT ii-c: Convergence & Computational Time (Wei)
    prices_vs_N        = zeros(size(N_vec));
    computational_times = zeros(size(N_vec));

    fprintf('\n--- CONVERGENCE & PERFORMANCE ANALYSIS (WEI) ---\n');
    for i = 1:length(N_vec)
        cur_N = N_vec(i);
        cur_h = TTM / cur_N;
        t_nodes_curr   = 0 : cur_h : TTM;
        treeDates_curr = refDate + days(round(365 * t_nodes_curr));
        tic;

        B0_curr = getDiscountFactorByZeroRatesLinearInterp(refDate, treeDates_curr, dates_dt, discounts);
        [x_curr, zeta_curr, kd_r_curr, ku_r_curr, p_r_curr] = buildXtree(cur_N, cur_h, kappa, sigma_r, zeta0);
        df_HW_curr = fwdDiscounts_OU(x_curr, kappa, sigma_r, t_nodes_curr, B0_curr);

        [S_curr, jd_curr, ju_curr, p_hat_curr] = buildStree(cur_N, cur_h, kappa, sigma_S, sigma_r, ...
                                                             S0, X0, rho, I, zeta_curr, yearFraction);
        prices_vs_N(i)        = backwardInduction(S_curr, cur_N, K_atm, p_r_curr, p_hat_curr, kd_r_curr, jd_curr, df_HW_curr);
        computational_times(i) = toc;
        fprintf('N = %d \t| Price = %.6f \t| Time = %.3f sec\n', cur_N, prices_vs_N(i), computational_times(i));
    end

    figure('Name', 'Point ii.c: Precision vs Time');
    yyaxis left;  plot(N_vec, prices_vs_N, '-ok', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); ylabel('ATM Option Price (USD)');
    yyaxis right; plot(N_vec, computational_times, '-^g', 'LineWidth', 1.5, 'MarkerFaceColor', 'g'); ylabel('Computational Time (s)');
    xlabel('Number of Time Steps (N)'); grid on;
    title('Convergence & Computational Cost vs N');

else  % method == 2: Robust

    %% POINT ii-a: American Put Price vs Strike (Robust)
    [S_tree_rob, q_ju_ku, q_ju_kd, q_jd_ku, q_jd_kd, ju_S_rob, jd_S_rob, ku_r_rob, kd_r_rob, df_rob] = ...
        tree_Robust(N, h, kappa, zeta0, sigma_S, sigma_r, S0, rho);

    prices_vs_K_rob = zeros(size(K_vec));
    for idx = 1:length(K_vec)
        prices_vs_K_rob(idx) = Robust_backward_Induction(S_tree_rob, N, K_vec(idx), ...
            q_ju_ku, q_ju_kd, q_jd_ku, q_jd_kd, ju_S_rob, jd_S_rob, ku_r_rob, kd_r_rob, df_rob);
                                                                       % ^^^^^^ was df_HW
    end

    figure('Name', 'Robust: American Put vs Strike');
    plot(K_vec, prices_vs_K_rob, '-o', 'LineWidth', 1.5); grid on;
    title('Robust - American Put Option Price vs Strike (T=2y, \rho=0.3)');
    xlabel('Strike Price (K)'); ylabel('Option Price (USD)');
    xline(S0, '--r', 'ATM Strike');
    legend('Robust American Put', 'ATM Level', 'Location', 'best');

    %% POINT ii-b: ATM Price vs Correlation (Robust)
     for i = 1:length(rho_new)
        current_rho = rho_new(i);
        [S_rob_r, q_uu_r, q_ud_r, q_du_r, q_dd_r, ju_r, jd_r, ku_rr, kd_rr, df_rob_r] = ...
            tree_Robust(N, h, kappa, zeta0, sigma_S, sigma_r, S0, current_rho, I, yearFraction);
        prices_vs_rho_rob(i) = Robust_backward_Induction(S_rob_r, N, K_atm, ...
            q_uu_r, q_ud_r, q_du_r, q_dd_r, ju_r, jd_r, ku_rr, kd_rr, df_rob_r);
    end

    figure('Name', 'Robust: American Put vs Correlation');
    plot(rho_new, prices_vs_rho_rob, '-ro', 'LineWidth', 1.5); grid on;
    title('Robust - ATM American Put Price vs Correlation \rho (T=2y, K=S_0)');
    xlabel('Correlation (\rho)'); ylabel('Option Price (USD)');

    %% POINT ii-c: Convergence & Computational Time (Robust)
    prices_vs_N_rob = zeros(size(N_vec));
    times_rob       = zeros(size(N_vec));

    fprintf('\n--- CONVERGENCE & PERFORMANCE ANALYSIS (ROBUST) ---\n');

    for i = 1:length(N_vec)
         cur_N = N_vec(i);
        cur_h = TTM / cur_N;
        t_nodes_curr   = 0 : cur_h : TTM;
        treeDates_curr = refDate + days(round(365 * t_nodes_curr));
        tic;
        [S_rob_c, q_uu_c, q_ud_c, q_du_c, q_dd_c, ju_c, jd_c, ku_c, kd_c, df_rob_c] = ...
            tree_Robust(cur_N, cur_h, kappa, zeta0, sigma_S, sigma_r, S0, rho, I, yearFraction);
        prices_vs_N_rob(i) = Robust_backward_Induction(S_rob_c, cur_N, K_atm, ...
            q_uu_c, q_ud_c, q_du_c, q_dd_c, ju_c, jd_c, ku_c, kd_c, df_rob_c);
        times_rob(i) = toc;
        fprintf('N = %d \t| Robust Price = %.6f \t| Time = %.3f sec\n', cur_N, prices_vs_N_rob(i), times_rob(i));
    end

    figure('Name', 'Robust: Precision vs Time');
    yyaxis left;  plot(N_vec, prices_vs_N_rob, '-ok', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); ylabel('ATM Option Price (USD)');
    yyaxis right; plot(N_vec, times_rob, '-^g', 'LineWidth', 1.5, 'MarkerFaceColor', 'g'); ylabel('Computational Time (s)');
    xlabel('Number of Time Steps (N)'); grid on;
    title('Robust - Convergence & Computational Cost vs N');


end



%% POINT iii): Volatility Calibration & True American Pricing
[S_tree, jd_S, ju_S, p_hat_S] = buildStree(N, h, kappa, sigma_S, sigma_r, ...
                                            S0, X0, rho, I, zeta_tree, yearFraction);
% 1. Retrieve specific T=2y target data from arrays
F_2y = fwd_mid_curve(9);      % 1.17480
df_2y = discounts(9);         % 0.97399953
vol_mkt = 0.1048;             % 10.48% (Market implied ATM Volatility)
K_atm = S0;                   % ATM-Spot

fprintf('\n--- POINT iii: VOLATILITY CALIBRATION & REPRICING ---\n');
fprintf('Target Market Volatility (Black): %.4f%%\n', vol_mkt * 100);

% Calculate the Target European Price for benchmarking
r_cont = -log(df_2y) / TTM;
[~, europeanPrice_Target] = blkprice(F_2y, K_atm, r_cont, TTM, vol_mkt);

% 2. Calibrate the true underlying stochastic volatility (sigma_S)
% using fzero to match the analytical European Stoch IR formula.
sigma_S_calibrated = calibrateVolatility(F_2y, TTM, vol_mkt, rho, kappa, sigma_r, K_atm, df_2y);
fprintf('Calibrated True Volatility (sigma_S): %.4f%%\n', sigma_S_calibrated * 100);

% 3. Reprice the American Option using the newly calibrated volatility
fprintf('\nRepricing American Option with calibrated volatility (N = %d)...\n', N);

% Rebuild FX Tree using sigma_S_calibrated instead of the raw market vol
[S_tree_calib, jd_S_calib, ju_S_calib, p_hat_S_calib] = buildStree(N, h, kappa, sigma_S_calibrated, sigma_r, ...
                                                                   S0, X0, rho, I, zeta_tree, yearFraction);

% Execute backward induction for the American option on the new tree
americanPrice_calib = backwardInduction(S_tree_calib, N, K_atm, p_r, p_hat_S_calib, kd_r, jd_S_calib, df_HW);

fprintf('American Option Price (CRR Tree):      %.6f\n', americanPrice_calib);
