% Mohsen Seifi, CS676, Winter 2021, University of Waterloo
% -------------------------------------------------------------------------
% risk-neutral pricing of down-and-out put option (closed form solution)
% v = cfdao(Price, Strike, Rate, Time, Maturity, Volatility, Barrier)
% Price: spot price of the underlying
% Strike: strike price of the down-and-out put option
% Rate: annualized interest rate 
% Time: price time
% Maturity: maturity of the down-and-out put option
% Volatility: annulaized percentage volatility of the underlying price
% Barrier: barrier of the down-and-out put option

function v = cfdao(Price, Strike, Rate, Time, Maturity, Volatility, ...
    Barrier)

sVect = Price;

k = Strike;
r = Rate;
t = Time;   T = Maturity;
sigma = Volatility; variance = sigma^2;

b = Barrier;

d = exp(-r * (T - t));

sk = log(sVect/k);
d1 = (sk + (r + variance/2) * sqrt(T-t))/(sigma);
d2 = (sk + (r - variance/2) * sqrt(T-t))/(sigma);

sb = log(sVect/b);
d3 = (sb + (r + variance/2) * sqrt(T-t))/(sigma);
d4 = (sb + (r - variance/2) * sqrt(T-t))/(sigma);
d5 = (sb - (r - variance/2) * sqrt(T-t))/(sigma);
d6 = (sb - (r + variance/2) * sqrt(T-t))/(sigma);

skb = log(sVect * k/b^2);
d7 = (skb - (r - variance/2) * sqrt(T-t))/(sigma);
d8 = (skb - (r + variance/2) * sqrt(T-t))/(sigma);

v = k * d * (normcdf(d4) - normcdf(d2) - (b./sVect).^ ...
    (-1 + 2 * r/variance) .* (normcdf(d7) - normcdf(d5))) - sVect .* ...
    (normcdf(d3) - normcdf(d1) - (b./sVect).^(1 + 2 * r/variance) .* ...
    (normcdf(d8) - normcdf(d6)));
end