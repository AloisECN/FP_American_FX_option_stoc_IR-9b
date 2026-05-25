function sigma = calibrateVolatility(F, yearFrac, sigma_S, rho, kappa, sigma_r, K, df)

T = yearFrac;
r = -1/yearFrac * log(df);
[~, putBlack]  = blkprice(F, K, r, T, sigma_S);

objective = @(sigma_opt) EuroPriceStochIR(F, T, sigma_opt, rho, kappa, sigma_r, K, df) - putBlack;

vol_mkt = sigma_S;
sigma = fzero(objective, vol_mkt);

end