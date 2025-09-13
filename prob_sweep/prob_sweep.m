% This matlab function is used to plot degradation-inducing effect of
% rejuvenation

clear all; close all;
params = initializeParams;
params.Nx = 400;
params.x = linspace(0,2,params.Nx+1);
params.dx = params.x(2) - params.x(1);

params.Tfinal = 400;
params.dt = 0.025;
params.Nt = params.Tfinal/params.dt;

params.vP = 0.05;
params.vW = 0.05;
params.lambdaR = 0;

params.delta = 0.6 * params.x';
params.k3 = 0;
params.k4 = 0;

lambdaPs = [1.08496,1.08443, 1.08404, 1.08376,1.08357];
p3s = [0.1,0.3,0.5,0.7,0.9];
k1s = [0.0001676, 0.00032969, 0.000505075,0.00077682,0.00158114]/params.dx;

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
    params.p3 = p3s(i);
    parama.p1 = (1-params.p3)/2;
    params.p2 = (1-params.p3)/2;

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

legend(h, arrayfun(@(r) sprintf('$p_3=%.1f$', r), p3s,... 
    'UniformOutput', false), 'Interpreter', 'latex', 'Location','best');
    title('Ratio Dynamics');

saveas(f1, 'ratio_p3', 'fig');
saveas(f1, 'ratio_p3', 'svg');

f2 = figure(2);
colors = lines(length(lambdaPs));
linestyles = {'-', '--', '-', '--', '-'};
h = gobjects(length(lambdaPs),1);

for i = 1:length(lambdaPs)
    params.lambdaP = lambdaPs(i);
    params.k1 = k1s(i);
    params.k2 = 0.1 * params.k1;
    params.p3 = p3s(i);
    parama.p1 = (1-params.p3)/2;
    params.p2 = (1-params.p3)/2;

    [PSol,WSol] = main(params);
    time_vec = (1:params.Nt+1)*params.dt;
    Pend = PSol(:,end);

    h(i) = plot(params.x(1:3*params.Nx/4+1), Pend(1:3*params.Nx/4+1), 'Color', colors(i,:), 'LineStyle', linestyles{i}); hold on;
end

xlabel('x', 'Interpreter', 'latex');
ylabel('Density', 'Interpreter', 'latex');
legend(h, arrayfun(@(r) sprintf('$p_3=%.1f$', r), p3s, ...
    'UniformOutput', false), 'Interpreter', 'latex', 'Location', 'best');
title('Damage Distribution for Stem', 'Interpreter','latex');

saveas(f2, 'damage_S_p3', 'fig');
saveas(f2, 'damage_S_p3', 'svg');

f3 = figure(3);
colors = lines(length(lambdaPs));
linestyles = {'-', '--', '-', '--', '-'};
h= gobjects(length(lambdaPs), 1);

for i = 1:length(lambdaPs)
    params.lambdaP = lambdaPs(i);
    params.k1 = k1s(i);
    params.k2 = 0.1 * params.k1;
    params.p3 = p3s(i);
    parama.p1 = (1-params.p3)/2;
    params.p2 = (1-params.p3)/2;

    [PSol,WSol] = main(params);
    time_vec = (1:params.Nt+1)*params.dt;
    Pend = PSol(:,end);
    Wend = WSol(:,end);

    h(i) = plot(params.x(1:3*params.Nx/4+1), Wend(1:3*params.Nx/4+1), 'Color', colors(i,:), 'LineStyle', linestyles{i}); hold on;
end

xlabel('x', 'Interpreter', 'latex');
ylabel('Density', 'Interpreter', 'latex');
legend(h, arrayfun(@(r) sprintf('$p_3=%.1f$', r), p3s, ...
    'UniformOutput', false), 'Interpreter', 'latex', 'Location','best');
title('Damage Distribution for TD', 'Interpreter','latex');

saveas(f3, 'damage_T_p3', 'fig');
saveas(f3, 'damage_T_p3', 'svg');