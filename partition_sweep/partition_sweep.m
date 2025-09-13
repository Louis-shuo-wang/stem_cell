% This matlab function is used to plot degradation-inducing effect of
% rejuvenation

clear all; close all;
params = initializeParams;
params.Nx = 400;
params.x = linspace(0,2,params.Nx+1);
params.dx = params.x(2) - params.x(1);

params.Tfinal = 2000;
params.dt = 0.025;
params.Nt = params.Tfinal/params.dt;

params.vP = 0.02;
params.vW = 0.02;
params.lambdaR = 0;

params.delta = 0.6 * params.x';
params.k3 = 0;
params.k4 = 0;

alphas = [0.1,0.3,0.5,0.5,0.5];
betas = [0.5,0.5,0.5,0.3,0.1];

lambdaPs = [0.6816,0.6837, 0.6842, 0.6825,0.6785];
k1s = [0.000009558, 0.0000095471,0.0000095389,0.0000095442,0.0000095545]/params.dx;

P_degra_ave = zeros(length(lambdaPs), 1);
W_degra_ave = zeros(length(lambdaPs), 1);
P_degra_max = zeros(length(lambdaPs), 1);
W_degra_max = zeros(length(lambdaPs), 1);
P_degra_right = zeros(length(lambdaPs), 1);
W_degra_right = zeros(length(lambdaPs), 1);
ratio_frac = zeros(length(lambdaPs), 1);

set(groot, 'defaultAxesFontSize', 30);
set(groot, 'defaultTextFontSize', 30);
set(groot, 'defaultLineLineWidth', 2.5);
set(groot, 'defaultAxesLineWidth', 1.2);
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesBox', 'off');

f1 = figure(1);
colors = lines(length(lambdaPs));
linestyles = {'-', '--', '-', '--', '-'};
h = gobjects(length(lambdaPs), 1);

for i = 1:length(lambdaPs)
    params.lambdaP = lambdaPs(i);
    params.k1 = k1s(i);
    params.k2 = 0.1 * params.k1;
    params.alpha1 = alphas(i);
    params.alpha2 = 1-params.alpha1;
    params.beta1 = betas(i);
    params.beta2 = 1-params.beta1;

    [PSol,WSol] = main(params);
    time_vec = (1:params.Nt+1)*params.dt;
    Psum = sum(PSol, 1) * params.dx;
    Wsum = sum(WSol, 1) * params.dx;
    ratio = Wsum ./ Psum;
    ratio_frac(i) = max(max(ratio - 7,0)/7);

    P_degra_sum = sum(linspace(0,2,params.Nx+1)'.*PSol, 1) * params.dx;
    P_degra_aves = P_degra_sum ./ Psum;
    W_degra_sum = sum(linspace(0,2,params.Nx+1)'.*WSol, 1) * params.dx;
    W_degra_aves = W_degra_sum ./ Wsum;
    P_degra_ave(i) = P_degra_aves(end);
    W_degra_ave(i) = W_degra_aves(end);

    [~,P_ind] = max(PSol(:, end));
    [~,W_ind] = max(WSol(:, end));
    P_degra_max(i) = P_ind * params.dx;
    W_degra_max(i) = W_ind * params.dx;
    P_right = find(PSol(:,end)>1e-4, 1, 'last')*params.dx;
    W_right = find(WSol(:,end)>1e-4, 1, 'last')*params.dx;
    P_degra_right(i) = P_right;
    W_degra_right(i) = W_right;

    h(i) = plot(time_vec(1:params.Nt/4+1), ratio(1:params.Nt/4+1), 'Color', colors(i,:), 'LineStyle', linestyles{i}); hold on;
end
 xlabel('t', 'Interpreter', 'latex');
    ylabel('Ratio', 'Interpreter', 'latex');

legend(h, arrayfun(@(ii) sprintf('$\\alpha_1=%.1f,\\,\\beta_1=%.1f$', ...
    alphas(ii), betas(ii)), 1:length(lambdaPs), 'UniformOutput', false), ...
    'Interpreter', 'latex', 'Location','best');

    title('Ratio Dynamics');

saveas(f1, 'ratio_alphabeta', 'fig');
saveas(f1, 'ratio_alphabeta', 'svg');

f2 = figure(2);
colors = lines(length(lambdaPs));
linestyles = {'-', '--', '-', '--', '-'};
h = gobjects(length(lambdaPs),1);

for i = 1:length(lambdaPs)
     params.lambdaP = lambdaPs(i);
    params.k1 = k1s(i);
    params.k2 = 0.1 * params.k1;
    params.alpha1 = alphas(i);
    params.alpha2 = 1-params.alpha1;
    params.beta1 = betas(i);
    params.beta2 = 1-params.beta1;

    [PSol,WSol] = main(params);
    time_vec = (1:params.Nt+1)*params.dt;
    Pend = PSol(:,end);

    h(i) = plot(params.x(1:3*params.Nx/4+1), Pend(1:3*params.Nx/4+1), 'Color', colors(i,:), 'LineStyle', linestyles{i}); hold on;
end

xlabel('x', 'Interpreter', 'latex');
ylabel('Density', 'Interpreter', 'latex');
legend(h, arrayfun(@(ii) sprintf('$\\alpha_1=%.1f,\\,\\beta_1=%.1f$', ...
    alphas(ii), betas(ii)), 1:length(lambdaPs), 'UniformOutput', false), ...
    'Interpreter', 'latex', 'Location','best');
title('Damage Distribution for Stem', 'Interpreter','latex');

saveas(f2, 'damage_S_alphabeta', 'fig');
saveas(f2, 'damage_S_alphabeta', 'svg');

f3 = figure(3);
colors = lines(length(lambdaPs));
linestyles = {'-', '--', '-', '--', '-'};
h= gobjects(length(lambdaPs), 1);

for i = 1:length(lambdaPs)
     params.lambdaP = lambdaPs(i);
    params.k1 = k1s(i);
    params.k2 = 0.1 * params.k1;
    params.alpha1 = alphas(i);
    params.alpha2 = 1-params.alpha1;
    params.beta1 = betas(i);
    params.beta2 = 1-params.beta1;

    [PSol,WSol] = main(params);
    time_vec = (1:params.Nt+1)*params.dt;
    Pend = PSol(:,end);
    Wend = WSol(:,end);

    h(i) = plot(params.x(1:3*params.Nx/4+1), Wend(1:3*params.Nx/4+1), 'Color', colors(i,:), 'LineStyle', linestyles{i}); hold on;
end

xlabel('x', 'Interpreter', 'latex');
ylabel('Density', 'Interpreter', 'latex');
legend(h, arrayfun(@(ii) sprintf('$\\alpha_1=%.1f,\\,\\beta_1=%.1f$', ...
    alphas(ii), betas(ii)), 1:length(lambdaPs), 'UniformOutput', false), ...
    'Interpreter', 'latex', 'Location','best');
title('Damage Distribution for TD', 'Interpreter','latex');

saveas(f3, 'damage_T_alphabeta', 'fig');
saveas(f3, 'damage_Tx_alphabeta', 'svg');