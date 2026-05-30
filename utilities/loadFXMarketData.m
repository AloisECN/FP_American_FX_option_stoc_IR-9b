function data = loadFXMarketData(filename)
% LOADFXMARKETDATA  Reads FX market data from the Excel file FXVolEURUSD.xlsm
%
%   DATA = LOADFXMARKETDATA(FILENAME)
%
%   Returned fields:
%       data.dates          - datetime (11x1) : USD discount curve dates
%       data.discounts      - double   (11x1) : USD discount factors
%       data.refDate        - datetime scalar : first pillar (spot date)
%       data.maturities     - datetime (10x1) : forward maturities (Spot->3Y)
%       data.fwd_mid_curve  - double   (10x1) : mid forward prices
%
%   Requires MATLAB R2019a or later.

    if nargin < 1
        filename = 'FXVolEURUSD.xlsm';
    end
    sheet = 'Values';
    
    % -- USD Discount Curve (G15:H25) 
    T_disc         = readtable(filename, 'Sheet', sheet, 'Range', 'G15:H25', ...
                               'ReadVariableNames', false);
    data.dates     = T_disc{:, 1};
    data.discounts = T_disc{:, 2};
    data.refDate   = data.dates(1);
    
    % -- Mid Forward Curve (G3:J12) --
    T_fwd              = readtable(filename, 'Sheet', sheet, 'Range', 'G3:J12', ...
                                   'ReadVariableNames', false);
    data.maturities    = T_fwd{:, 1};
    data.fwd_mid_curve = T_fwd{:, 4};
end