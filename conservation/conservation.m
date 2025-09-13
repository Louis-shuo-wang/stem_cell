% function that plot Figure 3;
clear all;close all;clf;
params = initializeParams;

params.Nx = 400;
params.x = linspace(0,4,params.Nx+1);
params.dx = params.x(2)-params.x(1);

params.lambdaP=1;
params.lambdaR=0;
params.delta = 0;

%%%%
params.p1=0.5; 
params.p2=0.5; 
params.p3=0;  
params.k1 = 0;
params.k2 = 0;
params.k3 = 0;

params.Tfinal = 20;
params.dt = 0.025;
params.Nt = params.Tfinal/params.dt;

params.vP = 0.025;
params.vW = 0.025;


[PSol, WSol] = main(params);
time_vec = (1:params.Nt+1)*params.dt;
Psum = sum(PSol,1)*params.dx;
Wsum = sum(WSol,1)*params.dx;

% damPsum = sum(PSol.*params.x')*params.dx;
% damWsum = sum(WSol.*params.x')*params.dx;
% damPave = damPsum ./ Psum;
% damWave = damWsum ./ Wsum;

set(groot, 'defaultAxesFontSize', 30);
set(groot, 'defaultTextFontSize', 30);
set(groot, 'defaultAxesLineWidth', 1.2);
set(groot, 'defaultLineLineWidth', 2.5);
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesBox', 'off');

figure(1);
plot(time_vec, Psum, 'b-', time_vec, Wsum, 'k--');
legend('$\hat{P}$','$\hat{W}$', 'Interpreter', 'latex', 'location', 'best');
xlabel('t', 'Interpreter','latex'); 
ylabel('$\mathrm{Total\,\,number}$', 'Interpreter', 'latex');
title('$v_P=v_W=0.025$', 'Interpreter','latex');


%%%%
% figure(1);
% plot(time_vec, Psum,'-r',time_vec, Wsum,'-b','LineWidth',2);
% legend('$\hat{P}$','$\hat{W}$', 'Interpreter', 'latex', 'Fontsize', 18, 'location', 'best');
% xlabel('Time', 'FontSize', 18); ylabel('Total number of units', 'FontSize', 18);
% title('Total number', 'FontSize', 18);
% 
% plot(params.x, PSol(:,end), 'b-', params.x, WSol(:,end), 'k--', 'LineWidth',2);
% legend('$\hat{P}$','$\hat{W}$', 'Interpreter', 'latex', 'Fontsize', 18, 'location', 'best');
% xlabel('x', 'FontSize', 18);
% ylabel('Density', 'FontSize', 18);
% title('Stationary Degradation Distribution', 'FontSize', 24);
% set(gca,'FontSize',18);
% % % After loop: show space-time images
% figure;
% for t = 1:10:params.Nt+1
% subplot(1,2,1);
% plot(params.x, PSol(:,t));
% axis xy; xlabel('x'); ylabel('P'); title(sprintf('P(a,%.2f)', (t-1)*params.dt));
% 
% subplot(1,2,2);
% plot(params.x, WSol(:,t));
% axis xy; xlabel('x'); ylabel('W'); title(sprintf('W(a,%.2f)', (t-1)*params.dt));
% pause(0.05);
% end


