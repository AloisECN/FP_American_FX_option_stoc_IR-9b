function [costOfCarry, I] = calibrateCostOfCarry(yearFrac, forwards, S0, sigma_r, k)
<<<<<<< HEAD
    
=======
    % CALIBRATECOSTOFCARRY Calibrates the cost-of-carry curve for an FX model
    % with stochastic interest rates (Extended Vasicek / Hull-White framework).
    %
    % INPUTS:
    %   yearFrac : Numeric vector. Market term structure (maturities T_k) in years.
    %   forwards : Numeric vector. Market Forward FX rates for each maturity.
    %   S0       : Scalar. Spot FX rate at t=0.
    %   sigma_r  : Scalar. Volatility of the domestic interest rate process.
    %   k        : Scalar. Mean reversion speed (kappa) of the domestic rate.
    %
    % OUTPUTS:
    %   costOfCarry : Numeric vector. Piecewise constant marginal rates extracted
    %                 via bootstrapping for each time interval.
    %   I           : Numeric vector. Total integral of the cost-of-carry from t=0
    %                 up to each market maturity T_k. Used for tree interpolation.

    % 1. Compute the analytical HJM convexity adjustment integral
    % We use element-wise operations (.^) to process all maturities simultaneously
>>>>>>> a365668 (Add MATLAB functions)
    E = exp(-k * yearFrac); 
    
    term1 = yearFrac;
    term2 = (2 / k) * (1 - E);
    term3 = (1 / (2 * k)) * (1 - E.^2); 
    
    hjm_int = (sigma_r^2 / k^2) * (term1 - term2 + term3);
    
<<<<<<< HEAD
    I = log(forwards ./ S0) + 0.5 * hjm_int;
    
    costOfCarry = zeros(size(yearFrac));
    
    costOfCarry(1) = I(1) / yearFrac(1);
    
=======
    % 2. Calculate the total integral I(T_k) of the cost-of-carry
    % Formula: I(T_k) = ln(F(0,T_k) / S0) + 0.5 * hjm_int
    I = log(forwards ./ S0) + 0.5 * hjm_int;
    
    % 3. Vectorized Bootstrapping to extract marginal cost-of-carry rates
    costOfCarry = zeros(size(yearFrac));
    
    % First bucket (from t=0 to T_1)
    costOfCarry(1) = I(1) / yearFrac(1);
    
    % Subsequent buckets (from T_{k-1} to T_k) 
    % diff(I) computes I(k) - I(k-1), diff(yearFrac) computes dt
>>>>>>> a365668 (Add MATLAB functions)
    costOfCarry(2:end) = diff(I) ./ diff(yearFrac);
    
end