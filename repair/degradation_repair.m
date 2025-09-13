% This function is used to plot degradation repair of rejuvenation

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
params.rho = 0.5;
params.delta = 0.6 * params.x';
params.p3 = 0;
params.k3 = 0;
params.k4 = 0;

lambdaPs = [1.08694, 1.08961, 1.09158, 1.09365, 1.09781];
lambdaRs = [0.01, 0.03, 0.05, 0.07, 0.09];
p1s = [0.5, 0.5, 0.5, 0.45, 0.3];
k1s = [0.000193532, 0.00040136, 0.00069716, 0.00110683, 0.0006339]/params.dx;

P_degra_ave = zeros(length(lambdaRs), 1);
W_degra_ave = zeros(length(lambdaRs), 1);
P_degra_max = zeros(length(lambdaRs), 1);
W_degra_max = zeros(length(lambdaRs), 1);
P_degra_right = zeros(length(lambdaRs), 1);
W_degra_right = zeros(length(lambdaRs), 1);

set(groot, 'defaultAxesFontSize', 30);
set(groot, 'defaultTextFontSize', 30);
set(groot, 'defaultLineLineWidth', 2.5);
set(groot, 'defaultAxesLineWidth', 1.2);
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesBox', 'off');

f1 = figure(1);
colors = lines(length(lambdaRs));
linestyles = {'-', '--', '-', '--', '-'};
h =gobjects(length(lambdaRs), 1);

for i = 1:length(lambdaRs)
    params.lambdaP = lambdaPs(i);
    params.lambdaR = lambdaRs(i);
    params.p1 = p1s(i);
    params.p2 = 1-params.p1;
    params.k1 = k1s(i);
    params.k2 = 0.1 * params.k1;

    [PSol, WSol] = main(params);
    time_vec = (1:params.Nt+1)*params.dt;
    Psum = sum(PSol, 1) * params.dx;
    Wsum = sum(WSol, 1) * params.dx;
    ratio = Wsum ./ Psum;
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
    P_right = find(PSol(:, end) > 1e-4, 1, 'last') * params.dx;
    W_right = find(WSol(:, end) > 1e-4, 1, 'last') * params.dx;
    P_degra_right(i) = P_right;
    W_degra_right(i) = W_right;

    h(i) = plot(time_vec, ratio, 'Color', colors(i, :), 'LineStyle', linestyles{i}); hold on;
end
xlabel('t', 'Interpreter', 'latex');
ylabel('Ratio', 'Interpreter', 'latex');

legend(h, arrayfun(@(r) sprintf('$\\lambda_R=%.2f$', r), lambdaRs, ...
    'UniformOutput', false), 'Interpreter', 'latex', 'Location', 'best');
title('Ratio Dynamics', 'Interpreter','latex');

saveas(f1, 'ratio_repair', 'fig');
saveas(f1, 'ratio_repair', 'svg');

f2 = figure(2);
colors = lines(length(lambdaRs));
linestyles = {'-', '--', '-', '--', '-'};
h = gobjects(length(lambdaRs), 1);

for i = 1:length(lambdaRs)
    params.lambdaP = lambdaPs(i);
    params.lambdaR = lambdaRs(i);
    params.p1 = p1s(i);
    params.p2 = 1 - params.p1;
    params.k1 = k1s(i);
    params.k2 = 0.1 * params.k1;

    [PSol, WSol] = main(params);
    time_vec = (1:params.Nt+1) * params.dt;
    Pend = PSol(:, end);

    h(i) = plot(params.x(1:3*params.Nx/4+1), Pend(1:3*params.Nx/4+1), 'Color', colors(i, :), 'LineStyle', linestyles{i}); hold on;
end

xlabel('x', 'Interpreter', 'latex');
ylabel('Density', 'Interpreter', 'latex');
legend(h, arrayfun(@(r) sprintf('$\\lambda_R=%.2f$', r), lambdaRs, ...
    'UniformOutput', false), 'Interpreter', 'latex', 'Location', 'best');
title('Damage Distribution for Stem', 'Interpreter', 'latex');

saveas(f2, 'damage_S_repair', 'fig');
saveas(f2, 'damage_S_repair', 'svg');

f3 = figure(3);
colors = lines(length(lambdaRs));
linestyles = {'-', '--', '-', '--', '-'};
h = gobjects(length(lambdaRs), 1);

for i = 1:length(lambdaRs)
    params.lambdaP = lambdaPs(i);
    params.lambdaR = lambdaRs(i);
    params.p1 = p1s(i);
    params.p2 = 1 - params.p1;
    params.k1 = k1s(i);
    params.k2 = 0.1 * params.k1;

    [PSol, WSol] = main(params);
    time_vec = (1:params.Nt+1) * params.dt;
    Wend = WSol(:, end);

    h(i) = plot(params.x(1:3*params.Nx/4+1), Wend(1:3*params.Nx/4+1), 'Color', colors(i, :), 'LineStyle', linestyles{i}); hold on;
end

xlabel('x', 'Interpreter', 'latex');
ylabel('Density', 'Interpreter', 'latex');
legend(h, arrayfun(@(r) sprintf('$\\lambda_R=%.2f$', r), lambdaRs, ...
    'UniformOutput', false), 'Interpreter', 'latex', 'Location', 'best');
title('Damage Distribution for TD', 'Interpreter','latex');

saveas(f3, 'damage_T_rapair', 'fig');
saveas(f3, 'damage_T_reapir', 'svg');
