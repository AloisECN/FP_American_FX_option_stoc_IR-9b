function sigma = calibrateVolatility(F, yearFrac, sigma_S, rho, kappa, sigma_r, K, df)
% CALIBRATEVOLATILITY Calibrates the spot volatility to match a market target price.
%
% This function finds the implied volatility parameter (sigma) for the 
% stochastic interest rate model that equates its European Put price 
% to the benchmark Black market price.
%
% INPUTS:
%   F        - Forward price of the underlying asset at maturity
%   yearFrac - Time to maturity (T) in years
%   sigma_S  - Market volatility (used to compute the target benchmark and as initial guess)
%   rho      - Correlation between the spot asset and the interest rate
%   kappa    - Mean reversion speed of the Hull-White interest rate model
%   sigma_r  - Volatility of the Hull-White interest rate model
%   K        - Strike price of the option
%   df       - Domestic discount factor from time 0 to maturity T
%
% OUTPUT:
%   sigma    - Calibrated spot volatility that zeroes out the pricing error

    T = yearFrac;
    
    % Convert the discount factor into an equivalent continuously compounded rate
    r = -1/yearFrac * log(df);
    
    %  2. Compute Benchmark Market Price 
    [~, putBlack]  = blkprice(F, K, r, T, sigma_S);
    
    objective = @(sigma_opt) EuroPriceStochIR(F, T, sigma_opt, rho, kappa, sigma_r, K, df) - putBlack;
    
    sigma = fzero(objective, sigma_S);
    
end