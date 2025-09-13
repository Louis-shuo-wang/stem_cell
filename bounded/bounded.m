% this is the function that plots Figure 3;

%  the left panel

clear all;close all;clf;

params = initializeParams;

params.Nx = 2000;
params.x = linspace(0,20,params.Nx+1);
params.dx = params.x(2)-params.x(1);

params.lambdaP=1;
params.lambdaR=0;
params.delta = 0.5;
% params.delta=0.6 * linspace(0,4,params.Nx+1)';

params.p1=0.5; 
params.p2=0.5; 
params.p3=0;  % base case
params.k1 = 0;
params.k2 = 0;
params.k3 = 0;

params.Tfinal = 200;
params.dt = 0.025;
params.Nt = params.Tfinal/params.dt;

params.vP = 0.2;
params.vW = 0.2;


[PSol, WSol] = main(params);
time_vec = (1:params.Nt+1)*params.dt;
Psum = sum(PSol,1)*params.dx;
Wsum = sum(WSol,1)*params.dx;

damPsum = sum(PSol.*params.x')*params.dx;
damWsum = sum(WSol.*params.x')*params.dx;
damPave = damPsum ./ Psum;
damWave = damWsum ./ Wsum;

set(groot, 'defaultAxesFontSize', 30);
set(groot, 'defaultTextFontSize', 30);
set(groot, 'defaultAxesLineWidth', 1.2);
set(groot, 'defaultLineLineWidth', 2.5);
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesBox', 'off');

% figure;
% plot(time_vec,damPave, 'b-', time_vec, damWave, 'k--');
% legend('P_ave_dam','W-ave_dam');


figure(1);
% plot(time_vec, Psum,'-r',time_vec, Wsum,'-b');
% legend('$\hat{P}$','$\hat{W}$', 'Interpreter', 'latex', 'location', 'best');
% xlabel('Time', 'Interpreter','latex'); ylabel('Total number of units');
% title('Total number', 'Interpreter','latex');

plot(params.x(1:params.Nx/5), PSol(1:params.Nx/5,end), 'b-', params.x(1:params.Nx/5), WSol(1:params.Nx/5,end));
legend('$P$','$W$', 'Interpreter', 'latex', 'location', 'best');
xlabel('x', 'Interpreter', 'latex');
ylabel('Density', 'Interpreter', 'latex');
title('Stationary Damage Distribution', 'Interpreter', 'latex');

saveas(gcf, 'bounded_a.fig');
saveas(gcf, 'bounded_b.svg')
% 
% figure;
% for t = 1:100:params.Nt+1
% subplot(1,2,1);
% plot(params.x, PSol(:,t));
% axis xy; xlabel('x'); ylabel('P'); title(sprintf('P(a,%.2f)', (t-1)*params.dt));
% 
% subplot(1,2,2);
% plot(params.x, WSol(:,t));
% axis xy; xlabel('x'); ylabel('W'); title(sprintf('W(a,%.2f)', (t-1)*params.dt));
% pause(0.05);
% end

