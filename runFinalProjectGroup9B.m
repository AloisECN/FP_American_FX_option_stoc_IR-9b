%% MAIN SCRIPT: Project - American FX Option with Stochastic Interest Rate
clc; clear all; close all;

%% 0. MARKET DATA INGESTION & INITIAL CALIBRATION

% Base Parameters
sigma_r = 0.01;
kappa = 0.10;
TTM = 2;
S0 = 1.14055;
sigma_S = 0.1048; % ATM Black volatility at T=2y
rho = 0.3;   % Base correlation

% Market Dates and Discount Factors (Domestic)
dates = { '12/Apr/2016'; '13/Apr/2016'; '19/Apr/2016'; '12/May/2016'; ...
          '13/Jun/2016'; '12/Jul/2016'; '12/Oct/2016'; '12/Apr/2017'; ...
          '12/Apr/2018'; '12/Apr/2019'; '14/Apr/2020'};

dates_dt = datetime(dates, 'InputFormat','dd/MMM/yyyy', 'Locale','en_US');
refDate = dates_dt(1);

discounts  = [ 1; 0.99998157; 0.99987514; 0.99942882; 0.99872386; ...
               0.99787021; 0.99480441; 0.98710348; 0.97399953; ...
               0.95687968; 0.93750149]; 

% Forward FX Mid Curve
fwd_mid_curve = [ 1.14055; 1.14077; 1.14152; 1.14257; 1.14370; ...
                  1.14730; 1.15135; 1.15558; 1.17480; 1.19685 ];

maturities = { '12/Apr/2016'; '19/Apr/2016'; '12/May/2016'; '13/Jun/2016'; ...
               '12/Jul/2016'; '12/Oct/2016'; '12/Jan/2017'; '12/Apr/2017'; ...
               '12/Apr/2018'; '12/Apr/2019'};

maturities_dt = datetime(maturities, 'InputFormat','dd/MMM/yyyy', 'Locale','en_US');

% Cost of Carry Calibration
yearFraction = yearfrac(refDate, maturities_dt, 3);
[costOfCarry, I] = calibrateCostOfCarry(yearFraction, fwd_mid_curve, S0, sigma_r, kappa);

% Base OU Process starting point (zeta0)
r0 = -1/(1/365)*log(0.99998157);
zeta0 = r0/sigma_r;

%% POINT ii) - a: American Put Prices vs Strike

N = 50; 
h = TTM / N;
t_nodes = 0 : h : TTM;
treeDates = refDate + days(round(365 * t_nodes));

% 1. Build Trees 
B0_nodes = getDiscountFactorByZeroRatesLinearInterp(refDate, treeDates, dates_dt, discounts);
[x_tree, zeta_tree, kd_r, ku_r, p_r] = buildXtree(N, h, kappa, sigma_r, zeta0);
df_HW = fwdDiscounts_OU(x_tree, kappa, sigma_r, t_nodes, B0_nodes);

X0 = -rho * zeta0 / (sqrt(1 - rho^2));
[S_tree, jd_S, ju_S, p_hat_S] = buildStree(N, h, kappa, sigma_S, sigma_r, ...
                                            S0, X0, rho, I, zeta_tree, yearFraction);

% 2. Loop over Strikes
K_vec = linspace(0.8 * S0, 1.2 * S0, 100); % 20 points for a smooth plot
prices_vs_K = zeros(size(K_vec));

for idx = 1:length(K_vec)
    prices_vs_K(idx) = backwardInduction(S_tree, N, K_vec(idx), p_r, p_hat_S, kd_r, jd_S, df_HW);
end

% 3. Plotting
figure('Name', 'American Put vs Strike');
plot(K_vec, prices_vs_K, '-ob', 'LineWidth', 1.5);
grid on;
title('American Put Option Price vs Strike (T=2y, \rho=0.3)');
xlabel('Strike Price (K)');
ylabel('Option Price (USD)');
xline(S0, '--r', 'ATM Strike');
legend('American Put Price', 'ATM Level', 'Location', 'best');


%% POINT ii) - b: ATM American Put Prices vs Correlation (rho)

K_atm = S0;
% Avoid exact -1 and 1 to prevent division by zero in X0 definition
rho_new = linspace(-0.99, 0.99, 50); 
prices_vs_rho = zeros(size(rho));

for i = 1:length(rho_new)
    current_rho = rho_new(i);
    
    % Recompute X0 and the S_tree (since the FX dynamics depend on rho)
    X0_current = -current_rho * zeta0 / sqrt(1 - current_rho^2);
    [S_tree_current, jd_S_curr, ju_S_curr, p_hat_S_curr] = buildStree(N, h, kappa, sigma_S, sigma_r, ...
                                                                     S0, X0_current, current_rho, I, zeta_tree, yearFraction);
    
    % Price the option (x_tree and df_HW do not change with rho)
    prices_vs_rho(i) = backwardInduction(S_tree_current, N, K_atm, p_r, p_hat_S_curr, kd_r, jd_S_curr, df_HW);
end

% Plotting
plot(rho_new, prices_vs_rho, '-or', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
grid on;
title('ATM American Put Price vs Correlation \rho (T=2y, K=S_0)');
xlabel('Correlation (\rho)');
ylabel('Option Price (USD)');


%% POINT ii) - c: Precision issues vs Computational Time

N_vec = [10, 25, 50, 75, 100, 150, 200, 300, 400];
prices_vs_N = zeros(size(N_vec));
computational_times = zeros(size(N_vec));

for i = 1:length(N_vec)
    current_N = N_vec(i);
    current_h = TTM / current_N;
    t_nodes_curr = 0 : current_h : TTM;
    treeDates_curr = refDate + days(round(365 * t_nodes_curr));
    
    tic; % Start timer
    
    % 1. Rebuild IR Tree for current N
    B0_nodes_curr = getDiscountFactorByZeroRatesLinearInterp(refDate, treeDates_curr, dates_dt, discounts);
    [x_tree_curr, zeta_tree_curr, kd_r_curr, ku_r_curr, p_r_curr] = buildXtree(current_N, current_h, kappa, sigma_r, zeta0);
    df_HW_curr = fwdDiscounts_OU(x_tree_curr, kappa, sigma_r, t_nodes_curr, B0_nodes_curr);
    
    % 2. Rebuild FX Tree for current N
    [S_tree_curr, jd_S_curr, ju_S_curr, p_hat_S_curr] = buildStree(current_N, current_h, kappa, sigma_S, sigma_r, ...
                                                                   S0, X0, rho, I, zeta_tree_curr, yearFraction);
                                                               
    % 3. Backward Induction
    prices_vs_N(i) = backwardInduction(S_tree_curr, current_N, K_atm, p_r_curr, p_hat_S_curr, kd_r_curr, jd_S_curr, df_HW_curr);
    
    computational_times(i) = toc; % Stop timer
    
    fprintf('N = %d | Price = %.6f | Time = %.3f sec\n', current_N, prices_vs_N(i), computational_times(i));
end

% Plotting dual-axis graph (Precision vs Time)
figure('Name', 'Point i.c: Precision vs Time');
yyaxis left
plot(N_vec, prices_vs_N, '-ok', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
ylabel('ATM Option Price (USD)');
yyaxis right
plot(N_vec, computational_times, '-^g', 'LineWidth', 1.5, 'MarkerFaceColor', 'g');
ylabel('Computational Time (seconds)');
xlabel('Number of Time Steps (N)');
grid on;
title('Convergence & Computational Cost vs N');

disp('--- Script Execution Completed ---');
