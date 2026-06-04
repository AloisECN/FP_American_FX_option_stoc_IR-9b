function sigma = calibrateVolatility(F, yearFrac, vol_mkt, rho, kappa, sigma_r, K, df)
% CALIBRATEVOLATILITY Calibrates sigma_S under stochastic IR
% so that the analytical European put matches the market Black price.

    T = yearFrac;

    % Market Black target price
    putBlack = BlackFXPutForward(F, K, df, T, vol_mkt);

    % Stochastic IR model price minus target
    objective = @(sigma_opt) EuroPriceStochIR(F, T, sigma_opt, rho, kappa, sigma_r, K, df) - putBlack;

    % Solve for calibrated spot volatility
    sigma = fzero(objective, vol_mkt);
end