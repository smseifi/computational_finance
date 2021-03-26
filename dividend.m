% handling dividend payment when pricing with binomial trees
% vBar = dividend(v, s, d)
% input(s):
% v:    value of option @ t^{+}
% s:    spot price of the underlying asset
% d:    discrete dollar dividend
% output(s): 
% vBar:    value of ption @ t^{-}

function vBar = dividend(v, s, d)

sBar = max(s-d, min(s));
vBar = interp1(s, v, sBar);
end