%% 1. (a)
clc
clear all
format compact
% ----------------------------------------------------------------------- %

n = 12; dZero = 5e-2; dVect = dZero ./ 2.^(0:1:n);

type = 0;
s = 100;    k = 95; r = 5e-2;   T = 1;  sigma = 18e-2;
vVect = zeros(n + 1, 1);
for i = 1:n+1
    d = dVect(i);
    vVect(i, 1) = binomprice(d, type, s, k, r, T, sigma);
end

[call, put] = blsprice(100, 95, 0.05, 1, 0.18);

if type == 1
    exact = call;
else
    exact = put;
end

% ----------------------------------------------------------------------- %

change = diff(vVect);   ratio = change(1:end-1) ./ change(2:end);
change = [nan; change];   ratio = [nan; nan; ratio];

format long
tb1 = table(dVect', vVect, change, ratio, 'VariableNames', ...
    {'Delta t', 'Value', 'Change', 'Ratio'});
disp(tb1)

%% 1. (b)
clc
clear all

% ----------------------------------------------------------------------- %

s = 100;    k = 95; r = 5e-2;   T = 1;  sigma = 18e-2; dZero = 1;
d = 5e-3;
rhoVect = [0, 2e-2, 5e-2, 1e-1];    rl = length(rhoVect);

vMat = zeros(2, rl);
for j = 1:rl
    rho = rhoVect(j);
    vMat(1, j) = binomprice(d, 1, s, k, r, T, sigma, 1, rho);
    vMat(2, j) = binomprice(d, 0, s, k, r, T, sigma, 1, rho);
end

% ----------------------------------------------------------------------- %

tb1Hat = table(rhoVect', vMat(1, :)', vMat(2, :)', 'VariableNames', ...
    {'rho', 'call', 'put'});
disp(tb1Hat)

%% 2. (a)
clc
clear all 

% ----------------------------------------------------------------------- %

sVect = 100:2:120;  k = 100;    b = 90;
t = 0;  T = 1;
sigma = 18e-2;  r = 5e-2; 

vVect = cfdao(sVect, k, r, t, T, sigma, b);

% ----------------------------------------------------------------------- %

f1 = figure(1);
[default_position, default_color] = plset;
f1.InnerPosition = default_position;
f1.OuterPosition = default_position;

plot(sVect, vVect, 'd:', 'Color', ...
    default_color(1, :), 'LineWidth', 3)

ax = gca;   ax.FontSize = 15;   ax.LineWidth = 1;

ax.XGrid = 'on';    ax.YGrid = 'on';
ax.XMinorTick = 'on';   ax.XMinorGrid = 'on';
ax.YMinorTick = 'on';   ax.YMinorGrid = 'on';
ax.YLim = [15e-2, 23e-2]; 

%% 2. (b) iii
clc
clear all

% ----------------------------------------------------------------------- %

deltaVect = [2e-2, 1e-2, 5e-3, 25e-4]; dl = length(deltaVect);
mVect = 1e+3:4e+3:45e+3; ml = length(mVect);

s = 100;  k = 100;    b = 90; fZero = 0;
t = 0;  T = 1;
sigma = 18e-2;  r = 5e-2; 

vMat = zeros(ml, dl);
for j = 1:dl
    d = deltaVect(1, j);
    for i = 1:ml
        m = mVect(1, i);
        vMat(i, j) = mcdao(s, k, r, t, T, sigma, b, fZero, d, m);
    end
end

% ----------------------------------------------------------------------- %

f2 = figure(2);
[default_position, default_color] = plset;
f2.InnerPosition = default_position;
f2.OuterPosition = default_position;

cfdao = 0.194942412804864;
plot(mVect, cfdao * ones(1, ml), 'k:', 'LineWidth', 2, ...
    'DisplayName', 'closed-form');
hold on 
for i = 1:dl
plot(mVect, vMat(:, i), '-o', 'Color',default_color(i, :), ...
    'LineWidth', 1.5, 'DisplayName', ...
    strcat('$\Delta t = $', num2str(deltaVect(1, i))))
hold on
end
hold off 

grid on
grid minor
legend('Location', 'southeast', 'Interpreter', 'LaTex', 'LineWidth', ...
    1.25, 'FontSize', 10);

ax = gca;   ax.FontSize = 15;   ax.LineWidth = 1;

ax.XGrid = 'on';
ax.XMinorTick = 'off';   ax.XMinorGrid = 'on';

ax.Units = 'normalized';


%% 2. (b) iv
clc
clear all

% ----------------------------------------------------------------------- %

s = 100;  k = 100;    b = 90;
t = 0;  T = 1;
sigma = 18e-2;  r = 5e-2;

fZeroVect = [0, 10, 20]; fl = length(fZeroVect);
delta = 25e-4;  m = 1e+5;
vVect = zeros(fl, 1);
for i = 1:fl
    f = fZeroVect(1, i);
    vVect(i, 1) = mcdao(s, k, r, t, T, sigma, b, f, delta, m);
end

% ----------------------------------------------------------------------- %

disp('problem 2. (b) iv:')

tb2 = table(fZeroVect', vVect, 'VariableNames', ...
    {'rebate', 'price'});
disp(tb2)

%% Graduate Student Question

clear all

% ----------------------------------------------------------------------- %

deltaVect = [2e-2, 1e-2, 5e-3, 25e-4]; dl = length(deltaVect);
m = 1e+5;

s = 100;  k = 100;    b = 90; fZero = 0;
t = 0;  T = 1;
sigma = 18e-2;  r = 5e-2; 

vVect = zeros(dl, 1);
for i = 1:dl
    d = deltaVect(1, i);
    vVect(i, 1) = mod_mcdao(s, k, r, t, T, sigma, b, fZero, d, m);
end

% ----------------------------------------------------------------------- %

disp('Graduate Student Question:')
tb3 = table(deltaVect', vVect, 'VariableNames', ...
    {'Delta t', 'price'});
disp(tb3)

%% Graduate Student Question / Comparison 
clc
clear all

% ----------------------------------------------------------------------- %

deltaVect = [2e-2, 1e-2, 5e-3, 25e-4]; dl = length(deltaVect);
mVect = 1e+3:4e+3:45e+3; ml = length(mVect);

s = 100;  k = 100;    b = 90; fZero = 0;
t = 0;  T = 1;
sigma = 18e-2;  r = 5e-2; 
tic
vMat = zeros(ml, dl);
for j = 1:dl
    d = deltaVect(1, j);
    for i = 1:ml
        m = mVect(1, i);
        vMat(i, j) = mod_mcdao(s, k, r, t, T, sigma, b, fZero, d, m);
    end
end
toc
% ----------------------------------------------------------------------- %

f3 = figure(3);
[default_position, default_color] = plset;
f3.InnerPosition = default_position;
f3.OuterPosition = default_position;

cfdao = 0.194942412804864;
plot(mVect, cfdao * ones(1, ml), 'k:', 'LineWidth', 2, ...
    'DisplayName', 'closed-form');
hold on 
for i = 1:dl
plot(mVect, vMat(:, i), '-o', 'Color',default_color(i, :), ...
    'LineWidth', 1.5, 'DisplayName', ...
    strcat('$\Delta t = $', num2str(deltaVect(1, i))))
hold on
end
hold off 

grid on
grid minor
legend('Location', 'southeast', 'Interpreter', 'LaTex', 'LineWidth', ...
    1.25, 'FontSize', 10);

ax = gca;   ax.FontSize = 15;   ax.LineWidth = 1;

ax.XGrid = 'on';
ax.XMinorTick = 'off';   ax.XMinorGrid = 'on';

ax.Units = 'normalized';

% ----------------------------------------------------------------------- %