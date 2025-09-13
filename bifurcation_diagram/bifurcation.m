% plot bifurcation diagram of (f, lambdaR)
% the relation between f and lambdaR is 2*delta*f = delta - lambdaR
% so this matlab function is used to visualize this relationship
clear all;close all;
lambdaR = linspace(0,1,101);
y1 = arrayfun(@(x) 0.5-1.25*x, lambdaR);
y2 = arrayfun(@(x) 0.5-x, lambdaR);
y3 = arrayfun(@(x) 0.5-5*x/6, lambdaR);

f1 = figure(1);
set(groot, 'defaultAxesFontSize', 24);
set(groot, 'defaultTextFontSize', 24);
set(groot, 'defaultLineLineWidth', 2.5);
set(groot, 'defaultAxesLineWidth', 1.2);
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesBox', 'off');

plot(lambdaR, y1, 'b-', lambdaR, y2, 'k--', lambdaR, y3, 'r');
xlabel('$\hat{\lambda}_R$','Interpreter','latex');
ylabel('$\hat{f}$', 'Interpreter','latex');
ylim([0,1]);
legend('$\delta=0.4$', '$\delta=0.5$', '$\delta=0.6$', 'Interpreter', 'latex', 'Location', 'northwest');
title('Bifurcation Diagram in the $(\hat{\lambda}_R, \hat{f})$ Plane', 'Interpreter', 'latex', 'FontSize', 21.5);
dim1 = [.2, .2, .2, .1];
dim2 = [.6, .7, .2, .1];
annotation('textbox', dim1, 'String', 'Extinction');
annotation('textbox', dim2, 'String', 'Growth');

saveas(f1, 'Bifurcation', 'fig');
saveas(f1, 'Bifurcation', 'svg');