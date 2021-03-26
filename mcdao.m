% Monte Carlo risk-neutral pricing of down-and-out put option 
% v = mcdao(Price, Strike, Rate, Time, Maturity, Volatility, Barrier, ...
% Rebate, Time_step, Paths)
% Price: spot price of the underlying
% Strike: strike price of the down-and-out put option
% Rate: annualized interest rate 
% Time: price time
% Maturity: maturity of the down-and-out put option
% Volatility: annulaized percentage volatility of the underlying price
% Barrier: barrier of the down-and-out put option
% Rebate: rebate of the down-and-out put option
% Time_step: time discretization mesh size
% Paths: number of simulated paths 

function v = mcdao(Price, Strike, Rate, Time, Maturity, Volatility, ...
    Barrier, Rebate, Time_step, Paths)

T = Maturity;   t = Time;   delta = Time_step;
n = (T - t)/delta;
m = Paths;

s0 = Price;    k = Strike; b = Barrier; fZero = Rebate;
sigma = Volatility; r = Rate;

v = zeros(m, 1);
for i = 1:m
    phiVect = randn(n, 1);
    rhs = (r - sigma^2/2) * delta + sigma * phiVect * sqrt(delta);
    sPath = [s0; exp(cumsum(rhs)+log(s0))];
    v(i, 1) = payoff(sPath, k, b, fZero);
end
v = exp(-r) * sum(v)/m;
end

function v = payoff(s, k, b, fZero)
if any(s - b <= 0)
    v = fZero;
else
    v = max(k - s(end), 0);
end
end