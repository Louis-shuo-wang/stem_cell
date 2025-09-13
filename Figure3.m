% function that plot Figure 3;
clear all;close all;clf;
params = initializeParams;
params.lambdaP=1;
params.lambdaR=0;
params.delta=0;

params.p1=0.5; 
params.p2=0.5; 
params.p3=0;  % base case


params.Tfinal = 25;
params.Nt = params.Tfinal/params.dt;

params.vP = 0;
params.vW = 0;

[PSol, WSol] = main(params); % each column is a set of solution on one time point
Psum = sum(PSol,1)*params.dx;
Wsum = sum(WSol,1)*params.dx;
timevec = (1:params.Nt+1)*params.dt;

set(groot, 'defaultAxesFontSize', 18);
set(groot, 'defaultTextFontSize', 18);
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultAxesLineWidth', 1.2);
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesBox', 'off');

figure('Position', [100, 100, 1200, 600]);
subplot(1,2,1);
plot(timevec, Psum, 'b-', timevec, Wsum, 'k--','LineWidth',1.5);
legend('$\hat{P}$', '$\hat{W}$', 'Interpreter', 'latex', 'FontSize', 15,'location','best');
xlabel('t', 'Interpreter', 'latex');
ylabel('Total Number', 'Interpreter', 'latex');
title('$v_P=v_W=0$', 'Interpreter', 'latex', 'FontSize', 24);
set(gca, 'FontSize', 20);
params.vP = 0.025;
params.vW = 0.025;

[PSol,WSol] = main(params);
Psum = sum(PSol,1)*params.dx;
Wsum = sum(WSol,1)*params.dx;

subplot(1,2,2);
plot(timevec,Psum,'b-',timevec,Wsum,'k--','LineWidth',1.5);
xlabel('t', 'FontSize', 18);
ylabel('Total Number', 'FontSize', 18);
title('$v_P=v_W=0.025$', 'Interpreter', 'latex','FontSize', 24);
set(gca, 'FontSize', 20);


