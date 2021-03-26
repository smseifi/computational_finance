% European option pricing using binomial trees (*handles dividend*) 
% vVect = binomprice(Time_step, Type, Price, Strike, Rate, Maturity, ...
%    Volatility, Dividend_floor, Dividend_factor)
% input(s):
% Time_step: time discretization mesh size
% Type = 1/0 indicates a European call/put
% Price: spot price of the underlying
% Strike: strike price of the option
% Rate: annualized interest rate 
% Maturity: maturity of the down-and-out put option
% Volatility: annulaized percentage volatility of the underlying price
% Dividend_floor: floor for dividend
% Dividend_factor: ratio of the spot price by which the stock pays dividend
% output(s):
% v: fair value of option @ time t = 0
% * in what follows, we have assumed that the dividend payment time is
% equal to (Maturity/4) *

function v = binomprice(Time_step, Type, Price, Strike, Rate, Maturity, ...
    Volatility, Dividend_floor, Dividend_factor)

sZero = Price;
k = Strike;
r = Rate;
T = Maturity;
sigma = Volatility;

delta = Time_step;

n = ceil(T/delta);
nD = n/4;   
if abs(floor(nD) - nD) <= .5
    nD = floor(nD);
else
    nD = ceill(nD);
end
    
discount = exp(-r * delta);

u = exp(sigma * sqrt(delta) + (r - (sigma^2)/2) * delta);
d = exp(-sigma * sqrt(delta) + (r - (sigma^2)/2) * delta);
q = (exp(r * delta) - d)/(u - d);

sVect = (sZero * u.^(0:n) .* d.^(n:-1:0)).';

switch Type
    case 1
        v = max(sVect - k, 0);
    case 0
        v = max(k - sVect, 0);
    otherwise
        type_msg = 'Type should be 1 for Call and 0 for Put';
        error(type_msg);
end
switch nargin
    case 7
        for i = n:-1:1
            v = discount * ((1-q)*v(1:i) + q*v(2:i+1));
        end
    case 9
        dZero = Dividend_floor;
        rho = Dividend_factor;
        
        for i = n:-1:nD
            sVect = discount * ((1-q)*sVect(1:i) + q*sVect(2:i+1));
            v = discount * ((1-q)*v(1:i) + q*v(2:i+1));
        end
        d = max(rho * sVect, dZero);
        v = dividend(v, sVect, d);
        
        for i = nD-1:-1:1
            v = discount * ((1-q)*v(1:i) + q*v(2:i+1));
        end
    otherwise
        arg_msg = 'Unexpected number of input arguments';
        error(arg_msg);
end
end