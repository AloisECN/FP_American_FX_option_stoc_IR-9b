function price = EuroPriceStochIR(F, yearFrac, sigma_S, rho, kappa, sigma_r, K, df)
    
    T = yearFrac;
    
    exp_kT  = exp(-kappa * T);
    exp_2kT = exp(-2 * kappa * T);
    
    var_S  = sigma_S^2 * T;
    
    cov_Sr = -2 * rho * sigma_S * (sigma_r / kappa) * (T - (1 - exp_kT) / kappa);
    
    var_r  = (sigma_r / kappa)^2 * (T - (2 / kappa) * (1 - exp_kT) + (1 / (2 * kappa)) * (1 - exp_2kT));
    
    Var = var_S + cov_Sr + var_r;
    
    sqrt_Var = sqrt(Var); 
    
    d1 = log(F / K) / sqrt_Var + 0.5 * sqrt_Var;
    d2 = d1 - sqrt_Var;
    
    price = df * (K * normcdf(-d2) - F * normcdf(-d1));
    
end