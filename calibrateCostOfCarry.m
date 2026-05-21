function [costOfCarry, I] = calibrateCostOfCarry(yearFrac, forwards, S0, sigma_r, k)
    
    E = exp(-k * yearFrac); 
    
    term1 = yearFrac;
    term2 = (2 / k) * (1 - E);
    term3 = (1 / (2 * k)) * (1 - E.^2); 
    
    hjm_int = (sigma_r^2 / k^2) * (term1 - term2 + term3);
    
    I = log(forwards ./ S0) + 0.5 * hjm_int;
    
    costOfCarry = zeros(size(yearFrac));
    
    costOfCarry(1) = I(1) / yearFrac(1);
    
    costOfCarry(2:end) = diff(I) ./ diff(yearFrac);
    
end