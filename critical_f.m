% this function is used to plot the critical renewal parameter
% f_{\rm{crit}} =
% \frac{\delta_W-\lambda_R}{2\delta_W}+\frac{\delta_W+\lambda_R}{2\delta_W\lambda_P}\delta_P
% and the original critical parameter without introducing the stem death
% f_{\rm{crit}} =
% \frac{\delta_W-\lambda_R}{2\delta_W}

deltaP = linspace(0,1,201);
deltaW = 0.5;
lambdaR = 0.2;
lambdaP = 1;

y1 = (deltaW-lambdaR)/(2*deltaW) * ones(size(deltaP));
y2 = arrayfun(@(x) (deltaW-lambdaR)/(2*deltaW) + (deltaW+lambdaR)/(2*deltaW*lambdaP)*x, deltaP);

colors = lines(2);
linestyles = {'-', '--'};

set(groot, 'defaultAxesFontSize', 30);
set(groot, 'defaultTextFontSize', 30);
set(groot, 'defaultAxesLineWidth', 1.2);
set(groot, 'defaultLineLineWidth', 2.5);
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesBox', 'off');

figure;

plot(deltaP, y1, 'Color', colors(1,:), 'LineStyle', linestyles{1},'Displayname','original $\tilde{f}$');
hold on;
plot(deltaP, y2, 'Color', colors(2,:), 'LineStyle', linestyles{2}, 'Displayname', 'new $\tilde{f}$');
xlabel('$\delta_P$', 'Interpreter', 'latex');
ylabel('$f_{\rm{crit}}$', 'Interpreter','latex');
legend('Interpreter', 'latex','Location','best');
title('Comparison of two critical thresholds', 'Interpreter','latex');