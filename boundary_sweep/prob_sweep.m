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

params.p1 = 0.5;
params.p2 = 1 - params.p1;
params.p3 = 0;

params.delta = 0.6 * params.x';
params.k3 = 0;
params.k4 = 0;

lambdaPs = [0.759348,0.826909,0.8933,0.958721,0.991066];
k1s = [0.09184,0.07052,0.05436,0.04035,0.0335];
gs = [5,4,3,2,1.5];

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

    idx = i;
    params.lambdaP = lambdaPs(i);
    params.k1 = k1s(i);
    params.k2 = 0.1 * params.k1;


    [PSol,WSol] = main(params, idx);
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

legend(h, arrayfun(@(r) sprintf('$g=%.1f$', r), gs,... 
    'UniformOutput', false), 'Interpreter', 'latex', 'Location','best');
    title('Ratio Dynamics');

saveas(f1, 'ratio_BC', 'fig');
saveas(f1, 'ratio_BC', 'svg');

f2 = figure(2);
colors = lines(length(lambdaPs));
linestyles = {'-', '--', '-', '--', '-'};
h = gobjects(length(lambdaPs),1);

for i = 1:length(lambdaPs)

    idx = i;

    params.lambdaP = lambdaPs(i);
    params.k1 = k1s(i);
    params.k2 = 0.1 * params.k1;


    [PSol,WSol] = main(params, idx);
    time_vec = (1:params.Nt+1)*params.dt;
    Pend = PSol(:,end);

    h(i) = plot(params.x(1:3*params.Nx/4+1), Pend(1:3*params.Nx/4+1), 'Color', colors(i,:), 'LineStyle', linestyles{i}); hold on;
end

xlabel('x', 'Interpreter', 'latex');
ylabel('Density', 'Interpreter', 'latex');
legend(h, arrayfun(@(r) sprintf('$g=%.1f$', r), gs, ...
    'UniformOutput', false), 'Interpreter', 'latex', 'Location', 'best');
title('Damage Distribution for Stem', 'Interpreter','latex');

saveas(f2, 'damage_S_BC', 'fig');
saveas(f2, 'damage_S_BC', 'svg');

f3 = figure(3);
colors = lines(length(lambdaPs));
linestyles = {'-', '--', '-', '--', '-'};
h= gobjects(length(lambdaPs), 1);

for i = 1:length(lambdaPs)

    idx = i;

    params.lambdaP = lambdaPs(i);
    params.k1 = k1s(i);
    params.k2 = 0.1 * params.k1;


    [PSol,WSol] = main(params, i);
    time_vec = (1:params.Nt+1)*params.dt;
    Pend = PSol(:,end);
    Wend = WSol(:,end);

    h(i) = plot(params.x(1:3*params.Nx/4+1), Wend(1:3*params.Nx/4+1), 'Color', colors(i,:), 'LineStyle', linestyles{i}); hold on;
end

xlabel('x', 'Interpreter', 'latex');
ylabel('Density', 'Interpreter', 'latex');
legend(h, arrayfun(@(r) sprintf('$g=%.1f$', r), gs, ...
    'UniformOutput', false), 'Interpreter', 'latex', 'Location','best');
title('Damage Distribution for TD', 'Interpreter','latex');

saveas(f3, 'damage_T_BC', 'fig');
saveas(f3, 'damage_T_BC', 'svg');