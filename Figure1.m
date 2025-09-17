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
% You can tweak parameters in the PARAMS section to explore regimes.
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
%% ===== Helper funcs =====
TwoP1minus1 = @(W) (2*hatp1 - 1) ./ (1 + (k1*W).^m1); % 2p1-1(W)
p1 = @(W) 0.5*(1 + TwoP1minus1(W)); % p1(W)
lambdaP = @(W) hatlambdaP ./ (1 + (k2*W).^m2); % lambda_P(W)
% Derivatives for the Jacobian at the equilibrium
% d(2p1-1)/dW = - (2*phat1 - 1) * [ m1*k1^m1 * W^(m1-1) ] / (1+(k1 W)^m1)^2
% Rewritten to be numerically stable
D_TwoP1minus1 = @(W) - (2*hatp1 - 1) .* (m1 .* (k1.^m1) .* (max(W,0).^(m1-1))) ./ (1 + (k1*W).^m1).^2;
D_lambdaP = @(W) - hatlambdaP .* (m2 .* (k2.^m2) .* (max(W,0).^(m2-1))) ./ (1 + (k2*W).^m2).^2;
% RHS of the ODE system
defRHS = @(t,Y) [ (TwoP1minus1(Y(2)).*lambdaP(Y(2))).*Y(1) + lambdaR*Y(2);
(2 - 2*p1(Y(2))).*lambdaP(Y(2)).*Y(1) - (delta + lambdaR).*Y(2) ];
Wscan = linspace(0,10,2001);
F = (1-2*p1(Wscan))*delta - lambdaR;
idx = find(sign(F(1:end-1)) .* sign(F(2:end)) <=0, 1);
% if isempty(idx)
% warning("No equilibrium");
% Wstar = NaN; Pstar = NaN;
% else
% Wstar = fzero(@(W) (1-2*p1(W))*delta-lambdaR, [Wscan(idx), Wscan(idx+1)]);
% Pstar = (delta/lambdaP(Wstar))*Wstar;
% end
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

set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 20);
set(groot, 'defaultTextFontSize', 20);
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultAxesLineWidth', 1.2);
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesBox', 'off');

%% ===== Vector field and nullclines =====
figure('Color','w','Position',[60 60 900 700]); hold on; box on;


%% --- Vector field (稀疏网格，用于箭头) ---
Ng_quiver = 25;
Pg_q = linspace(0,Pmax,Ng_quiver);
Wg_q = linspace(0,Wmax,Ng_quiver);
[PPq, WWq] = meshgrid(Pg_q, Wg_q);
UUq = (TwoP1minus1(WWq).*lambdaP(WWq)).*PPq + lambdaR.*WWq; % dP/dt
VVq = (2 - 2*p1(WWq)).*lambdaP(WWq).*PPq - (delta + lambdaR).*WWq; % dW/dt

% Normalize arrows
mag = sqrt(UUq.^2 + VVq.^2);
max_mag = max(mag);
mag = floor(1000 * (mag ./ max_mag))/1000;
mag(mag==0) = 1;
UUq = UUq./mag;
VVq = VVq./mag;

h1 = quiver(PPq,WWq,UUq,VVq,0.8,'Color',[0.6 0.6 0.6]);

%% --- Nullclines (细网格，用于 contour) ---
Ng_contour = 1000;
Pg_c = linspace(0,Pmax,Ng_contour);
Wg_c = linspace(0,Wmax,Ng_contour);
[PPc, WWc] = meshgrid(Pg_c, Wg_c);
UUc = (TwoP1minus1(WWc).*lambdaP(WWc)).*PPc + lambdaR.*WWc; % dP/dt
VVc = (2 - 2*p1(WWc)).*lambdaP(WWc).*PPc - (delta + lambdaR).*WWc; % dW/dt

[~, h2] = contour(PPc,WWc,UUc,[0 0],'LineWidth',2,'LineColor',[0.1 0.4 0.95]);
[~, h3] = contour(PPc,WWc,VVc,[0 0],'LineWidth',2,'LineColor',[0.85 0.2 0.2]);


legend([h1 h2 h3], {'Vector field','$d\hat{P}/dt = 0$','$d\hat{W}/dt = 0$'},'Location','best');
% Labels
xlabel('$\hat P$', 'Interpreter', 'latex'); ylabel('$\hat W$', 'Interpreter', 'latex'); 
title('P/W Totals: Phase Portrait with Saddle & Separatrix','Interpreter','none');
axis([0 Pmax 0 Wmax]); 
%% ===== Plot equilibrium and separatrix =====
if ~isnan(Wstar)
h_saddle = plot(Pstar,Wstar,'kp','MarkerSize',12,'MarkerFaceColor','y','DisplayName', 'Saddle');
% Build tiny displacements along eigenvectors
eps_sep = 0.1; % small offset along eigenvectors
y0_stable_1 = [Pstar;Wstar] + eps_sep * (v_stable / norm(v_stable));
y0_stable_2 = [Pstar;Wstar] - eps_sep * (v_stable / norm(v_stable));
y0_unstable_1 = [Pstar;Wstar] + eps_sep * (v_unstable / norm(v_unstable));
y0_unstable_2 = [Pstar;Wstar] - eps_sep * (v_unstable / norm(v_unstable));
% Integrate forward/backward to trace stable/unstable manifolds
eventStop = @(t,y) stopOnBounds(t,y,Pmax,Wmax); % stop when leaving plot window
odeoptsSep = odeset(odeopts,'Events',eventStop);
[~,Yb1] = ode45(@(t,y) -defRHS(t,y), [0 Tmax], y0_stable_1, odeoptsSep);
[~,Yb2] = ode45(@(t,y) -defRHS(t,y), [0 Tmax], y0_stable_2, odeoptsSep);
h_stable = plot([Yb1(:,1); NaN; Yb2(:,1)], [Yb1(:,2); NaN; Yb2(:,2)], ...
'k-','DisplayName','Stable manifold');
[~,Yf1] = ode45(defRHS, [0 10*Tmax], y0_unstable_1, odeoptsSep);
[~,Yf2] = ode45(defRHS, [0 10*Tmax], y0_unstable_2, odeoptsSep);
h_unstable = plot([Yf1(:,1); NaN; Yf2(:,1)], [Yf1(:,2); NaN; Yf2(:,2)], ...
'k--','DisplayName','Unstable manifold');
text(Pstar, Wstar, ' saddle','FontWeight','bold');
end
%% ===== Sample trajectories to illustrate fates =====
% Initial conditions sprinkled around the plane
ICs = [ 1 0.3; 3, 0.5; 1, 4; 3,4];
colors = lines(size(ICs,1));
h_traj = gobjects(size(ICs,1),1);
for i=1:size(ICs,1)
y0 = ICs(i,:)';
[T,Y] = ode45(defRHS,[0 10*Tmax],y0,odeopts);
% Plot trajectory line
h_traj(i) = plot(Y(:,1),Y(:,2),'-','LineWidth', 1.5, 'Color',colors(i,:),'DisplayName', sprintf('Traj%d',i));
% % ===== Add arrows along trajectory =====
% step = round(length(T)/10); % 10个箭头左右，更少
% idx = 1:step:length(T)-1;
% U = Y(idx+1,1)-Y(idx,1); % ΔP
% V = Y(idx+1,2)-Y(idx,2); % ΔW
% scale = 0.8; % 箭头稍大
% quiver(Y(idx,1),Y(idx,2),U,V,scale, ...
% 'Color',colors(i,:), 'MaxHeadSize',1.2, ...
% 'AutoScale','on','AutoScaleFactor',1.2, ...
% 'HandleVisibility','off'); % 不进入图例
end

legend([h1 h2 h3 h_stable h_unstable h_traj(:)'], ...
    {'Vector field','$d\hat{P}/dt = 0$','$d\hat{W}/dt = 0$', ...
     'Stable manifold','Unstable manifold', ...
     'Traj1','Traj2','Traj3','Traj4'}, ...
    'Interpreter','latex','FontSize',18,'Location','best');

hold off;
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