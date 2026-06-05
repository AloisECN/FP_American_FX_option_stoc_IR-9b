%% MAIN SCRIPT: Project - American FX Option with Stochastic Interest Rate
clc; clear all; close all;

%% 0. MARKET DATA INGESTION & INITIAL CALIBRATION
d = loadFXMarketData('FXVolEURUSD.xlsm');

dates_dt      = d.dates;           
refDate       = d.refDate;         
discounts     = d.discounts;       
fwd_mid_curve = d.fwd_mid_curve;   
maturities_dt = d.maturities;

sigma_r  = 0.01;
kappa    = 0.10;
S0       = 1.14055;
sigma_S  = 0.1048;
rho      = 0.3;

maturity_target = maturities_dt(9);
TTM = yearfrac(refDate, maturity_target, 3);

yearFraction = yearfrac(refDate, maturities_dt, 3);
[costOfCarry, I] = calibrateCostOfCarry(yearFraction, fwd_mid_curve, S0, sigma_r, kappa);
plotCostOfCarry(yearFraction, costOfCarry);

r0    = -1/(1/365) * log(discounts(2));
zeta0 = r0 / sigma_r;

%% Shared parameters (used by both methods and by Point iii))
N     = 300;
h     = TTM / N;
t_nodes  = (0:N)*h;
treeDates = refDate + days(round(365 * t_nodes));

K_atm    = S0;
K_vec    = linspace(0.8 * S0, 1.2 * S0, 100);
rho_new  = linspace(-0.99, 0.99, 100);
N_vec    = [10, 25, 50, 100, 150, 200, 250, 300, 400, 500];

% IR tree and discounts — needed by both methods
B0_nodes = getDiscountFactorByZeroRatesLinearInterp(refDate, treeDates, dates_dt, discounts);
[x_tree, zeta_tree, kd_r, ku_r, p_r] = buildXtree(N, h, kappa, sigma_r, zeta0);
df_HW    = fwdDiscounts_OU(x_tree, kappa, sigma_r, t_nodes, B0_nodes);
Y0       = -rho * zeta0 / sqrt(1 - rho^2);

%% Method selection
method = input('Choose the method: 1 --> Wei  |  2 --> Robust: ');

switch method
    
    case 1 %  WEI / STANDARD METHOD 
        fprintf('\n --- EXECUTING WEI METHOD --- \n');
        
        %% POINT ii)-a: American Put Price vs Strike (Wei)
        % The tree structure does not depend on the strike price K.
        % We build the tree only once outside the loop to optimize performance.
        [S_tree, jd_S, ju_S, p_hat_S] = buildStree(N, h, kappa, sigma_S, sigma_r, ...
                                                    S0, Y0, rho, I, zeta_tree, yearFraction);
        prices_vs_K = zeros(size(K_vec));
        for idx = 1:length(K_vec)
            prices_vs_K(idx) = backwardInduction(S_tree, zeta_tree, N, K_vec(idx), p_r, p_hat_S, ...
                kd_r, ku_r, jd_S, ju_S, df_HW, rho, sigma_r, sigma_S, h, 'Wei');
        end
        
        plotPriceVsStrike(K_vec, prices_vs_K, S0);
        
        %% POINT ii)-b: ATM Price vs Correlation (Wei)
        prices_vs_rho = zeros(size(rho_new));
        
        % CRITICAL: In Wei's orthogonalization, the transformed starting node Y0 
        % explicitly depends on the correlation rho. Therefore, the spatial 
        % tree MUST be rebuilt entirely for every new value of rho.
        for i = 1:length(rho_new)
            current_rho   = rho_new(i);
            Y0_current    = -current_rho * zeta0 / sqrt(1 - current_rho^2);
            
            [S_tree_curr, jd_S_curr, ju_S_curr, p_hat_S_curr] = buildStree(N, ...
                h, kappa, sigma_S, sigma_r, S0, Y0_current, current_rho, I, zeta_tree, yearFraction);
            
            prices_vs_rho(i) = backwardInduction(S_tree_curr, zeta_tree, N, K_atm, p_r, p_hat_S_curr, ...
                kd_r, ku_r, jd_S_curr, ju_S_curr, df_HW, current_rho, sigma_r, sigma_S, h, 'Wei');
        end
        
        plotPriceVsCorrelation(rho_new, prices_vs_rho);
        
        %% POINT ii)-c: Convergence & Computational Time (Wei)
        prices_vs_N         = zeros(size(N_vec));
        computational_times = zeros(size(N_vec));
        
        fprintf('\n --- CONVERGENCE & PERFORMANCE ANALYSIS (WEI / STANDARD) --- \n');
        
        % As the number of steps N changes, the time step h changes as well.
        % We must completely reconstruct the timeline, the forward discounts, 
        % and both the interest rate and FX trees from scratch.
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
                                                                 S0, Y0, rho, I, zeta_curr, yearFraction);
            
            prices_vs_N(i) = backwardInduction(S_curr, zeta_curr, cur_N, K_atm, p_r_curr, p_hat_curr, ...
                                               kd_r_curr, ku_r_curr, jd_curr, ju_curr, df_HW_curr, ...
                                               rho, sigma_r, sigma_S, cur_h, 'Wei');
            
            computational_times(i) = toc;
            
            fprintf('N = %d \t| Price = %.6f \t| Time = %.3f sec\n', cur_N, prices_vs_N(i), computational_times(i));
        end
        
        plotPrecisionVsTime(N_vec, prices_vs_N, computational_times, 'Wei');
        plotErrorConvergence(N_vec, prices_vs_N, 'Wei');

    case 2 %  ROBUST METHOD (Appolloni et al.) 
        fprintf('\n --- EXECUTING ROBUST METHOD ---\n');
        
        %% POINT ii)-a: American Put Price vs Strike (Robust)
        % Build the robust tree once, independently of the strike price K.
        U0 = log(S0) / sigma_S;
        [S_tree_rob, jd_S_rob, ju_S_rob, p_hat_S_rob] = buildStreeRobust(N, h, sigma_S, U0, x_tree, I, yearFraction);
        %%
        prices_vs_K_rob = zeros(size(K_vec));
        for idx = 1:length(K_vec)
            prices_vs_K_rob(idx) = backwardInduction(S_tree_rob, x_tree, N, K_vec(idx), p_r, p_hat_S_rob, ...
                kd_r, ku_r, jd_S_rob, ju_S_rob, df_HW, rho, sigma_r, sigma_S, h, 'Robust');
        end
        
        plotPriceVsStrike(K_vec, prices_vs_K_rob, S0);
        
        %% POINT ii)-b: ATM Price vs Correlation (Robust)
        prices_vs_rho_rob = zeros(size(rho_new));
        
        % Unlike Wei, the Robust spatial grid is independent of rho, so we 
        % only update probabilities during backward induction without rebuilding the tree.
        for i = 1:length(rho_new)
            current_rho = rho_new(i);
            
            prices_vs_rho_rob(i) = backwardInduction(S_tree_rob, x_tree, N, K_atm, p_r, p_hat_S_rob, ...
                                             kd_r, ku_r, jd_S_rob, ju_S_rob, df_HW, ...
                                             current_rho, sigma_r, sigma_S, h, 'Robust');
        end
        
        plotPriceVsCorrelation(rho_new, prices_vs_rho_rob);
        
        %% POINT ii)-c: Convergence & Computational Time (Robust)
        prices_vs_N_rob = zeros(size(N_vec));
        times_rob       = zeros(size(N_vec));
        
        fprintf('\n CONVERGENCE & PERFORMANCE ANALYSIS (ROBUST) \n');
        
        % Rebuilding the full architecture for each different N step size.
        for i = 1:length(N_vec)
            cur_N = N_vec(i);
            cur_h = TTM / cur_N;
            t_nodes_curr   = 0 : cur_h : TTM;
            treeDates_curr = refDate + days(round(365 * t_nodes_curr));
            
            tic; 
            
            B0_curr = getDiscountFactorByZeroRatesLinearInterp(refDate, treeDates_curr, dates_dt, discounts);
            [x_curr, zeta_curr, kd_r_curr, ku_r_curr, p_r_curr] = buildXtree(cur_N, cur_h, kappa, sigma_r, zeta0);
            df_HW_curr = fwdDiscounts_OU(x_curr, kappa, sigma_r, t_nodes_curr, B0_curr);
            
            U0_curr = log(S0) / sigma_S;
            [S_rob_c, jd_c, ju_c, p_hat_c] = buildStreeRobust(cur_N, cur_h, sigma_S, U0_curr, x_curr, I, yearFraction);
            
            prices_vs_N_rob(i) = backwardInduction(S_rob_c, x_curr, cur_N, K_atm, p_r_curr, p_hat_c, ...
                                           kd_r_curr, ku_r_curr, jd_c, ju_c, df_HW_curr, ...
                                           rho, sigma_r, sigma_S, cur_h, 'Robust');
                
            times_rob(i) = toc; 
            
            fprintf('N = %d \t| Robust Price = %.6f \t| Time = %.3f sec\n', cur_N, prices_vs_N_rob(i), times_rob(i));
        end
        
        plotPrecisionVsTime(N_vec, prices_vs_N_rob, times_rob, 'Robust');
        plotErrorConvergence(N_vec, prices_vs_N_rob, 'Robust');
end

%% POINT iii): Volatility Calibration & True American Pricing

F_2y = fwd_mid_curve(9);      % 1.17480
df_2y = discounts(9);         % 0.97399953
vol_mkt = 0.1048;             % 10.48% 
K_atm = S0;                   % ATM-Spot

fprintf('\n POINT iii: VOLATILITY CALIBRATION & REPRICING \n');
fprintf('Target Market Volatility (Black): %.4f%%\n', vol_mkt * 100);

r_cont = -log(df_2y) / TTM;
[~, europeanPrice_Target] = blkprice(F_2y, K_atm, r_cont, TTM, vol_mkt);

sigma_S_calibrated = calibrateVolatility(F_2y, TTM, vol_mkt, rho, kappa, sigma_r, K_atm, df_2y);
fprintf('Calibrated True Volatility (sigma_S): %.4f%%\n', sigma_S_calibrated * 100);

fprintf('\nRepricing American Option with calibrated volatility (N = %d)...\n', N);

switch method
    case 1 % WEI Method
        [S_tree_calib, jd_S_calib, ju_S_calib, p_hat_S_calib] = buildStree(N, h, kappa, sigma_S_calibrated,...
            sigma_r, S0, Y0, rho, I, zeta_tree, yearFraction);
            
        americanPrice_calib = backwardInduction(S_tree_calib, x_tree, N, K_atm, p_r, p_hat_S_calib, ...
                                                kd_r, ku_r, jd_S_calib, ju_S_calib, df_HW, ...
                                                rho, sigma_r, sigma_S_calibrated, h, 'Wei');
        methodName = 'Wei';
        
    case 2 % ROBUST Method (Appolloni et al.)
        U0_calib = log(S0) / sigma_S_calibrated;
        
        [S_tree_calib, jd_S_calib, ju_S_calib, p_hat_S_calib] = buildStreeRobust(N, h, sigma_S_calibrated, U0_calib, x_tree, I, yearFraction);
        
        americanPrice_calib = backwardInduction(S_tree_calib, x_tree, N, K_atm, p_r, p_hat_S_calib, ...
                                                kd_r, ku_r, jd_S_calib, ju_S_calib, df_HW, ...
                                                rho, sigma_r, sigma_S_calibrated, h, 'Robust');
        methodName = 'Robust';
 end

fprintf('American Option Price (%s Tree): %.6f\n', methodName, americanPrice_calib);

%% POINT iii) -d: Cost-of-Carry Parallel Shifts Analysis
shift = [-0.02, -0.01, 0.01, 0.02];

price_am_shifted = zeros(size(shift));
price_eur_shifted = zeros(size(shift));

fprintf('\nShift (%%) | European | American\n');
fprintf('\n');

switch method
    case 1
        for i = 1 : length(shift)
            s = shift(i); 
            
            % The function internally adjusts the integral of the Cost-of-Carry (drift).
            [S_tree_shifted, jd_S_shifted, ju_S_shifted, p_hat_S_shifted] = buildStree(N, h, kappa, sigma_S_calibrated, ...
                           sigma_r, S0, Y0, rho, I, zeta_tree, yearFraction, s);
                
            % Price the American Put option using the newly shifted tree
            price_am_shifted(i) = backwardInduction(S_tree_shifted, x_tree, N, K_atm, p_r, p_hat_S_shifted, ...
                                                kd_r, ku_r, jd_S_shifted, ju_S_shifted, df_HW, ...
                                                rho, sigma_r, sigma_S_calibrated, h, 'Wei');
            
            
            %  2. EUROPEAN OPTION (Analytical Method) 
            % the target Forward rate must be shifted consistently with the Cost-of-Carry.
            F_shifted = F_2y * exp(s * TTM);
            
            % Compute the closed-form European Put price using the shifted Forward
            price_eur_shifted(i) = EuroPriceStochIR(F_shifted, TTM, sigma_S_calibrated,...
                                                    rho, kappa, sigma_r, K_atm, df_2y);
                                                    
            fprintf('Shift = %+3.0f%% \t| Eur Price = %.5f \t| Am Price = %.5f\n', ...
                s * 100, price_eur_shifted(i), price_am_shifted(i));
        end

    case 2
        for i = 1 : length(shift)
            s = shift(i); 
            
            % The function internally adjusts the integral of the Cost-of-Carry (drift).
            [S_tree_shifted, jd_S_shifted, ju_S_shifted, p_hat_S_shifted] = buildStreeRobust(N, h, sigma_S_calibrated, U0, x_tree, I, yearFraction, s);       
            % Price the American Put option using the newly shifted tree
            price_am_shifted(i) = backwardInduction(S_tree_shifted, x_tree, N, K_atm, p_r, p_hat_S_shifted, ...
                                                kd_r, ku_r, jd_S_shifted, ju_S_shifted, df_HW, ...
                                                rho, sigma_r, sigma_S_calibrated, h, 'Robust');
            
            
            %  2. EUROPEAN OPTION (Analytical Method) 
            % the target Forward rate must be shifted consistently with the Cost-of-Carry.
            F_shifted = F_2y * exp(s * TTM);
            
            % Compute the closed-form European Put price using the shifted Forward
            price_eur_shifted(i) = EuroPriceStochIR(F_shifted, TTM, sigma_S_calibrated,...
                                                    rho, kappa, sigma_r, K_atm, df_2y);
                                                    
            fprintf('Shift = %+3.0f%% \t| Eur Price = %.5f \t| Am Price = %.5f\n', ...
                s * 100, price_eur_shifted(i), price_am_shifted(i));
        end
end


