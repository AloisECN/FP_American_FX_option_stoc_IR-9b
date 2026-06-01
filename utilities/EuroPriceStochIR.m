function price = EuroPriceStochIR(F, T, sigma_S, rho, kappa, sigma_r, K, df)
% EUROPRICESTOCHIR Prices a European Put option with stochastic interest rates.
% 
% This function computes the closed-form price of a European Put option
% assuming the underlying asset follows a lognormal diffusion and the 
% domestic interest rate follows a Hull-White (Extended Vasicek) process.
%
% INPUTS:
%   F        - Forward price of the underlying asset at maturity T
%   T        - Time to maturity (T) in years
%   sigma_S  - Volatility of the underlying spot asset
%   rho      - Correlation between the spot asset and the interest rate
%   kappa    - Mean reversion speed of the Hull-White interest rate model
%   sigma_r  - Volatility of the Hull-White interest rate model
%   K        - Strike price of the option
%   df       - Domestic discount factor from time 0 to maturity T (B_d(0,T))
%
% OUTPUT:
%   price    - Present value of the European Put option
    % 1. Pre-compute exponential terms for the Hull-White integrals 
    exp_kT  = exp(-kappa * T);
    exp_2kT = exp(-2 * kappa * T);
    
    % 2. Compute Integrated Variance Components (integrated from 0 to T)
    var_S  = sigma_S^2 * T;
    cov_Sr = rho * sigma_S * (sigma_r / kappa) * (T - (1 - exp_kT) / kappa);
    var_r  = (sigma_r / kappa)^2 * (T - (2 / kappa) * (1 - exp_kT) + (1 / (2 * kappa)) * (1 - exp_2kT));
    
    % Total integrated variance \Sigma^2(0,T) under the T-forward measure
    Var = var_S - 2*cov_Sr + var_r;
    sqrt_Var = sqrt(Var); 
    
    % 3. Black-like Pricing Formula 
    % d1 and d2 correspond directly to the derived text values
    d1 = (log(F / K) + 0.5 * Var) / sqrt_Var;
    d2 = d1 - sqrt_Var;
    
    % Final Put option pricing execution
    price = df * (K * normcdf(-d2) - F * normcdf(-d1));
    
end