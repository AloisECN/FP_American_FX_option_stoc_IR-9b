function [price, d1, d2, Sigma, Var] = EuroPriceStochIR(F, T, sigma_S, rho, kappa, sigma_r, K, df)
% EUROPRICESTOCHIR Prices a European FX put under stochastic interest rates.
%
% Formula:
%   P = df * (K * N(-d2) - F * N(-d1))
%
% with stochastic-IR adjusted total variance:
%   Var = sigma_S^2*T - 2*rho*sigma_S*IntSigmaB + IntSigmaB2

    if F <= 0 || K <= 0 || df <= 0 || T <= 0
        error('Inputs F, K, df and T must be positive.');
    end

    exp_kT  = exp(-kappa * T);
    exp_2kT = exp(-2 * kappa * T);

    % FX variance
    var_S = sigma_S^2 * T;

    % Integral of bond volatility sigma(u,T)
    intSigmaB = (sigma_r / kappa) * ...
        (T - (1 - exp_kT) / kappa);

    % Integral of bond volatility squared
    intSigmaB2 = (sigma_r / kappa)^2 * ...
        (T - (2 / kappa) * (1 - exp_kT) ...
        + (1 / (2 * kappa)) * (1 - exp_2kT));

    % Total integrated variance
    Var = var_S - 2 * rho * sigma_S * intSigmaB + intSigmaB2;

    if Var < 0 && abs(Var) < 1e-12
        Var = 0;
    elseif Var < 0
        error('Negative total variance encountered: Var = %.12g', Var);
    end

    Sigma = sqrt(Var);

    if Sigma == 0
        price = df * max(K - F, 0);
        d1 = NaN;
        d2 = NaN;
        return;
    end

    d1 = (log(F / K) + 0.5 * Var) / Sigma;
    d2 = d1 - Sigma;

    price = df * (K * normcdf(-d2) - F * normcdf(-d1));
end