%% Progenitor–Worker (P/W) totals: phase portrait, saddle & separatrix
% This script reproduces the open-loop stability picture in the paper
% for the *totals* (\hat P, \hat W). We use constant \delta and \lambda_R,
% while p1(\hat W) and \lambda_P(\hat W) have negative feedback.
%
% dP/dt = (2p1-1)*lambdaP*P + lambdaR*W
% dW/dt = 2(1-p1)*lambdaP*P - (delta + lambdaR)*W
%
% The equilibrium (P*, W*) satisfies two relations:
% (1) (1 - 2 p1(W*)) * delta = lambdaR
% (2) P*/W* = delta / lambdaP(W*)
% We solve (1) for W* and then use (2) to get P*.
%
% The script plots:
% • Vector field (quiver)
% • Nullclines (dP=0 and dW=0)
% • Saddle point (if found) and its stable/unstable separatrix
% • Sample trajectories that either go to extinction or diverge
%
% Publication-quality enhancements:
% - Increased line widths and marker sizes
% - Improved color scheme for better differentiation
% - LaTeX formatting for all labels
% - High-resolution output settings
% - Clean, minimal aesthetic with proper axis labeling
% - Consistent style throughout

clear; clc; close all;

%% ===== PARAMS =====
% Failure and rejuvenation (constants)
delta = 0.2; % Worker failure (>=0)
lambdaR = 0.02; % Rejuvenation from W->P (>=0)

% Feedback on p1(W) via 2p1-1 = (2*phat1 - 1) / (1 + (k1 W)^m1)
% (Pick phat1 < 0.5 so that 2p1-1 < 0 at W=0.)
hatp1 = 0.25; % base prob of symmetric self-renewal (no feedback)
k1 = 1; % feedback strength on p1
m1 = 2; % Hill exponent on p1

% Feedback on lambdaP(W) = lambdahatP / (1 + (k2 W)^m2)
hatlambdaP = 1.00; % base replication rate (no feedback)
k2 = 1; % feedback strength on lambdaP
m2 = 2; % Hill exponent on lambdaP

% ODE integration options
Tmax = 200; % max time for trajectories
odeopts = odeset('RelTol',1e-8,'AbsTol',1e-10);
Pmax = 4;
Wmax = 4;

%% ===== Publication Quality Figure Settings =====
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaultLineLineWidth',1.5);
set(0,'defaultAxesFontSize',12);
set(0,'defaultTextFontSize',14);

% Color scheme for publication
colors = struct();
colors.nullcline1 = [0, 0.4470, 0.7410]; % Blue for dP/dt=0
colors.nullcline2 = [0.8500, 0.3250, 0.0980]; % Orange for dW/dt=0
colors.stable = [0.4660, 0.6740, 0.1880]; % Green for stable manifold
colors.unstable = [0.6350, 0.0780, 0.1840]; % Red for unstable manifold
colors.saddle = [0.4940, 0.1840, 0.5560]; % Purple for saddle point
colors.trajectory = [0.3010, 0.7450, 0.9330; % Light blue trajectories
                    0.9290, 0.6940, 0.1250; % Yellow
                    0.6350, 0.0780, 0.1840; % Red
                    0.4660, 0.6740, 0.1880]; % Green

%% ===== Helper funcs =====
TwoP1minus1 = @(W) (2*hatp1 - 1) ./ (1 + (k1*W).^m1); % 2p1-1(W)
p1 = @(W) 0.5*(1 + TwoP1minus1(W)); % p1(W)
lambdaP = @(W) hatlambdaP ./ (1 + (k2*W).^m2); % lambda_P(W)

% Derivatives for the Jacobian at the equilibrium
D_TwoP1minus1 = @(W) - (2*hatp1 - 1) .* (m1 .* (k1.^m1) .* (max(W,0).^(m1-1))) ./ (1 + (k1*W).^m1).^2;
D_lambdaP = @(W) - hatlambdaP .* (m2 .* (k2.^m2) .* (max(W,0).^(m2-1))) ./ (1 + (k2*W).^m2).^2;

% RHS of the ODE system
defRHS = @(t,Y) [ (TwoP1minus1(Y(2)).*lambdaP(Y(2))).*Y(1) + lambdaR*Y(2);
                 (2 - 2*p1(Y(2))).*lambdaP(Y(2)).*Y(1) - (delta + lambdaR).*Y(2) ];

Wscan = linspace(0,10,2001);
F = (1-2*p1(Wscan))*delta - lambdaR;
idx = find(sign(F(1:end-1)) .* sign(F(2:end)) <=0, 1);

Wstar = (((1-2*hatp1)*delta - lambdaR)/lambdaR)^(1/m1)/k1;
Pstar = delta * Wstar / lambdaP(Wstar);

%% ===== Jacobian and eigen-structure at equilibrium =====
if ~isnan(Wstar)
    a = TwoP1minus1(Wstar); % (2p1 - 1) at W*
    ap = D_TwoP1minus1(Wstar); % d(2p1-1)/dW at W*
    b = lambdaP(Wstar); % lambdaP at W*
    bp = D_lambdaP(Wstar); % d lambdaP / dW at W*

    J = [ a*b, lambdaR + bp*a*Pstar + b*ap*Pstar;
         (2-2*p1(Wstar))*b, -(delta+lambdaR) + bp*(2-2*p1(Wstar))*Pstar - b*ap*Pstar ];

    [V,D] = eig(J);
    eigvals = diag(D);
    [~,i_pos] = max(real(eigvals));
    [~,i_neg] = min(real(eigvals));
    v_unstable = V(:,i_pos);
    v_stable = V(:,i_neg);
else
    J = NaN(2); eigvals = [NaN;NaN]; v_unstable = [NaN;NaN]; v_stable = [NaN;NaN];
end

%% ===== Create Publication Quality Figure =====
fig = figure('Color','w','Position',[100 100 800 600]);
hold on; box on;
set(gca,'Layer','top'); % Ensure axes are on top of plots

%% ===== Vector field =====
Ng = 20; % Reduced grid points for cleaner appearance
Pg = linspace(0,Pmax,Ng);
Wg = linspace(0,Wmax,Ng);
[PP,WW] = meshgrid(Pg,Wg);

UU = (TwoP1minus1(WW).*lambdaP(WW)).*PP + lambdaR.*WW; % dP/dt
VV = (2 - 2*p1(WW)).*lambdaP(WW).*PP - (delta + lambdaR).*WW; % dW/dt

% Normalize arrows
mag = sqrt(UU.^2 + VV.^2);
max_mag = max(mag(:));
scale_factor = 0.8/max_mag; % Adjusted scale factor
UU = UU * scale_factor;
VV = VV * scale_factor;

% Plot vector field with improved styling
h_quiver = quiver(PP, WW, UU, VV, 0.6, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, ...
                 'MaxHeadSize', 0.8, 'AutoScale', 'off');

%% ===== Nullclines =====
% Create finer grid for nullclines
Ng_fine = 100;
Pg_fine = linspace(0, Pmax, Ng_fine);
Wg_fine = linspace(0, Wmax, Ng_fine);
[PP_fine, WW_fine] = meshgrid(Pg_fine, Wg_fine);

UU_fine = (TwoP1minus1(WW_fine).*lambdaP(WW_fine)).*PP_fine + lambdaR.*WW_fine;
VV_fine = (2 - 2*p1(WW_fine)).*lambdaP(WW_fine).*PP_fine - (delta + lambdaR).*WW_fine;

% Plot nullclines with thicker lines
[~, h_dP] = contour(PP_fine, WW_fine, UU_fine, [0 0], 'LineWidth', 2.5, ...
                   'LineColor', colors.nullcline1, 'LineStyle', '-');
[~, h_dW] = contour(PP_fine, WW_fine, VV_fine, [0 0], 'LineWidth', 2.5, ...
                   'LineColor', colors.nullcline2, 'LineStyle', '-');

%% ===== Plot equilibrium and separatrix =====
if ~isnan(Wstar)
    % Plot saddle point
    h_saddle = plot(Pstar, Wstar, 'o', 'MarkerSize', 10, 'MarkerFaceColor', colors.saddle, ...
                   'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'DisplayName', 'Saddle Point');

    % Build tiny displacements along eigenvectors
    eps_sep = 0.05;
    y0_stable_1 = [Pstar; Wstar] + eps_sep * (v_stable / norm(v_stable));
    y0_stable_2 = [Pstar; Wstar] - eps_sep * (v_stable / norm(v_stable));
    y0_unstable_1 = [Pstar; Wstar] + eps_sep * (v_unstable / norm(v_unstable));
    y0_unstable_2 = [Pstar; Wstar] - eps_sep * (v_unstable / norm(v_unstable));

    % Integrate to trace stable/unstable manifolds
    eventStop = @(t,y) stopOnBounds(t,y,Pmax,Wmax);
    odeoptsSep = odeset(odeopts, 'Events', eventStop);

    % Stable manifold (backward integration)
    [~, Yb1] = ode45(@(t,y) -defRHS(t,y), [0 Tmax], y0_stable_1, odeoptsSep);
    [~, Yb2] = ode45(@(t,y) -defRHS(t,y), [0 Tmax], y0_stable_2, odeoptsSep);
    h_stable = plot([Yb1(:,1); NaN; Yb2(:,1)], [Yb1(:,2); NaN; Yb2(:,2)], ...
                   '-', 'LineWidth', 2.5, 'Color', colors.stable, 'DisplayName', 'Stable Manifold');

    % Unstable manifold (forward integration)
    [~, Yf1] = ode45(defRHS, [0 10*Tmax], y0_unstable_1, odeoptsSep);
    [~, Yf2] = ode45(defRHS, [0 10*Tmax], y0_unstable_2, odeoptsSep);
    h_unstable = plot([Yf1(:,1); NaN; Yf2(:,1)], [Yf1(:,2); NaN; Yf2(:,2)], ...
                     '--', 'LineWidth', 2.5, 'Color', colors.unstable, 'DisplayName', 'Unstable Manifold');

    % Label saddle point
    text(Pstar + 0.1, Wstar + 0.1, 'Saddle', 'FontSize', 12, 'FontWeight', 'bold', ...
         'Interpreter', 'latex', 'Color', colors.saddle);
end

%% ===== Sample trajectories =====
ICs = [1 0.3; 3, 0.5; 1, 4; 3, 4];
h_traj = gobjects(size(ICs,1),1);

for i = 1:size(ICs,1)
    y0 = ICs(i,:)';
    [T,Y] = ode45(defRHS, [0 10*Tmax], y0, odeopts);

    % Plot trajectory with arrow markers
    h_traj(i) = plot(Y(:,1), Y(:,2), '-', 'LineWidth', 2, 'Color', colors.trajectory(i,:));

    % Add direction arrows at selected points
    arrow_indices = round(linspace(1, length(T)-1, 5));
    for j = arrow_indices
        arrow_start = Y(j,:);
        arrow_end = Y(j+1,:);
        arrow(arrow_start, arrow_end, 'Length', 8, 'BaseAngle', 60, ...
              'TipAngle', 30, 'Color', colors.trajectory(i,:), 'Width', 1.5);
    end
end

%% ===== Final formatting =====
xlabel('Progenitor Population $\hat{P}$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Worker Population $\hat{W}$', 'Interpreter', 'latex', 'FontSize', 14);
title('Phase Portrait: Progenitor--Worker System', 'Interpreter', 'latex', 'FontSize', 16);

axis([0 Pmax 0 Wmax]);
grid on;
set(gca, 'GridAlpha', 0.2, 'MinorGridAlpha', 0.1);

% Create legend
legend_handles = [h_dP, h_dW];
legend_labels = {'$\frac{d\hat{P}}{dt} = 0$', '$\frac{d\hat{W}}{dt} = 0$'};

if ~isnan(Wstar)
    legend_handles = [legend_handles, h_stable, h_unstable, h_saddle];
    legend_labels = [legend_labels, {'Stable Manifold', 'Unstable Manifold', 'Saddle Point'}];
end

legend(legend_handles, legend_labels, 'Interpreter', 'latex', ...
       'Location', 'best', 'FontSize', 11, 'Box', 'off');

% Set aspect ratio to equal for proper phase portrait
axis equal;

% Improve overall appearance
set(gca, 'TickDir', 'out', 'TickLength', [0.02 0.02], ...
         'XMinorTick', 'on', 'YMinorTick', 'on', ...
         'LineWidth', 1.2);

hold off;

%% ===== Print parameters =====
disp('--- PARAMETERS USED ---');
fprintf('delta = %.4g, lambdaR = %.4g\n', delta, lambdaR);
fprintf('phat1 = %.3f, k1 = %.4g, m1 = %g\n', hatp1, k1, m1);
fprintf('lambdahatP = %.4g, k2 = %.4g, m2 = %g\n', hatlambdaP, k2, m2);

if ~isnan(Wstar)
    fprintf('Equilibrium (saddle) at: P* = %.6g, W* = %.6g\n', Pstar, Wstar);
    fprintf('Eigenvalues(J): [% .4e, % .4e]\n', eigvals(1), eigvals(2));
else
    fprintf('No equilibrium found in the plotting box. Adjust parameters.\n');
end

%% ===== Event: stop integration when leaving plotting box =====
function [value,isterminal,direction] = stopOnBounds(~,y,Pmax,Wmax)
    P = y(1); W = y(2);
    out = (P < 0) || (P > Pmax) || (W < 0) || (W > Wmax);
    value = 1 - double(out); % 0 when out of bounds
    isterminal = 1; % stop the integration
    direction = 0;
end

%% ===== Custom arrow function for trajectory direction markers =====
function arrow(start, stop, varargin)
    % Parse input parameters
    p = inputParser;
    addParameter(p, 'Length', 10, @isnumeric);
    addParameter(p, 'BaseAngle', 60, @isnumeric);
    addParameter(p, 'TipAngle', 30, @isnumeric);
    addParameter(p, 'Color', [0 0 0], @isnumeric);
    addParameter(p, 'Width', 1, @isnumeric);
    parse(p, varargin{:});

    % Calculate arrow direction and length
    direction = stop - start;
    length_dir = norm(direction);
    if length_dir == 0
        return;
    end
    direction = direction / length_dir;

    % Calculate perpendicular direction
    perp = [-direction(2); direction(1)];

    % Calculate arrow points
    tip = stop;
    base1 = start + direction * (length_dir - p.Results.Length * cosd(p.Results.BaseAngle)) + ...
            perp * (p.Results.Length * sind(p.Results.BaseAngle));
    base2 = start + direction * (length_dir - p.Results.Length * cosd(p.Results.BaseAngle)) - ...
            perp * (p.Results.Length * sind(p.Results.BaseAngle));

    % Draw arrow
    patch([tip(1), base1(1), base2(1)], [tip(2), base1(2), base2(2)], ...
          p.Results.Color, 'EdgeColor', p.Results.Color, 'LineWidth', p.Results.Width);
end
