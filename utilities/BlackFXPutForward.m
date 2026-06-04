function [price, d1, d2, Sigma] = BlackFXPutForward(F, K, df, T, vol)
% BLACKFXPUTFORWARD Prices a European FX put using the forward Black formula.
%
% Formula:
%   P = df * (K * N(-d2) - F * N(-d1))
%
% where:
%   Sigma = vol * sqrt(T)
%   d1 = (log(F/K) + 0.5*Sigma^2) / Sigma
%   d2 = d1 - Sigma

    if F <= 0 || K <= 0 || df <= 0 || T <= 0 || vol < 0
        error('Inputs F, K, df, T must be positive and vol must be non-negative.');
    end

    Sigma = vol * sqrt(T);

    if Sigma == 0
        price = df * max(K - F, 0);
        d1 = NaN;
        d2 = NaN;
        return;
    end

    d1 = (log(F / K) + 0.5 * Sigma^2) / Sigma;
    d2 = d1 - Sigma;

    price = df * (K * normcdf(-d2) - F * normcdf(-d1));
end