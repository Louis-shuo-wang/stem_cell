% function that plot Figure 3;
clear all;close all;clf;
params = initializeParams;

params.Nx = 400;
params.x = linspace(0,4,params.Nx+1);
params.dx = params.x(2)-params.x(1);

params.lambdaP=1;
params.lambdaR=0.05;
params.delta = 0.1;
% params.delta=0.6 * linspace(0,4,params.Nx+1)';



set(groot, 'defaultAxesFontSize', 30);
set(groot, 'defaultTextFontSize', 30);
set(groot, 'defaultAxesLineWidth', 1.2);
set(groot, 'defaultLineLineWidth', 2.5);
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesBox', 'off');


%%%
clear all;

params = initializeParams;

params.Nx = 400;
params.x = linspace(0,4,params.Nx+1);
params.dx = params.x(2)-params.x(1);

params.lambdaP=1;
params.lambdaR=0.05;
params.delta = 0.1;

params.p1=0.3; 
params.p2=1-params.p1; 
params.p3=0;  
params.k1 = 0;
params.k2 = 0;
params.k3 = 0;

params.Tfinal = 200;
params.dt = 0.025;
params.Nt = params.Tfinal/params.dt;

params.vP = 0.05;
params.vW = 0.05;


[PSol, WSol] = main(params);
time_vec = (1:params.Nt+1)*params.dt;
Psum = sum(PSol,1)*params.dx;
Wsum = sum(WSol,1)*params.dx;

% damPsum = sum(PSol.*params.x')*params.dx;
% damWsum = sum(WSol.*params.x')*params.dx;
% damPave = damPsum ./ Psum;
% damWave = damWsum ./ Wsum;

f1 = figure(1);
plot(time_vec, log(Psum), 'b-', time_vec, log(Wsum), 'k--');
legend('$\hat{P}$','$\hat{W}$', 'Interpreter', 'latex', 'location', 'best');
xlabel('t', 'Interpreter','latex'); 
ylabel('$\log(\mathrm{Total\,\,number})$', 'Interpreter', 'latex');
title('$\hat{f}=0.30>\tilde{f}$', 'Interpreter','latex');

pP = polyfit(time_vec(2001:8001), log(Psum(2001:8001)), 1);
pW = polyfit(time_vec(2001:8001), log(Wsum(2001:8001)), 1);

saveas(f1, 'bifurcation_a','fig');
saveas(f1, 'bifurcation_a','svg');

%%%%
clear all;

params = initializeParams;

params.Nx = 400;
params.x = linspace(0,4,params.Nx+1);
params.dx = params.x(2)-params.x(1);

params.lambdaP=1;
params.lambdaR=0.05;
params.delta = 0.1;

params.p1=0.25; 
params.p2=1-params.p1; 
params.p3=0;  
params.k1 = 0;
params.k2 = 0;
params.k3 = 0;

params.Tfinal = 200;
params.dt = 0.025;
params.Nt = params.Tfinal/params.dt;

params.vP = 0.05;
params.vW = 0.05;


[PSol, WSol] = main(params);
time_vec = (1:params.Nt+1)*params.dt;
Psum = sum(PSol,1)*params.dx;
Wsum = sum(WSol,1)*params.dx;

% damPsum = sum(PSol.*params.x')*params.dx;
% damWsum = sum(WSol.*params.x')*params.dx;
% damPave = damPsum ./ Psum;
% damWave = damWsum ./ Wsum;

f2 = figure(2);
plot(time_vec, Psum, 'b-', time_vec, Wsum, 'k--');
legend('$\hat{P}$','$\hat{W}$', 'Interpreter', 'latex', 'location', 'best');
xlabel('t', 'Interpreter','latex'); 
ylabel('$\mathrm{Total\,\,number}$', 'Interpreter', 'latex');
title('$\hat{f}=0.25=\tilde{f}$', 'Interpreter','latex');

% pP = polyfit(time_vec(2001:8001), Psum(2001:8001), 1);
% pW = polyfit(time_vec(2001:8001), Wsum(2001:8001), 1);

saveas(f2, 'bifurcation_b', 'fig');
saveas(f2, 'bifurcation_b', 'svg');

%%%%
clear all;

params = initializeParams;

params.Nx = 400;
params.x = linspace(0,4,params.Nx+1);
params.dx = params.x(2)-params.x(1);

params.lambdaP=1;
params.lambdaR=0.05;
params.delta = 0.1;

params.p1=0.2; 
params.p2=1-params.p1; 
params.p3=0;  
params.k1 = 0;
params.k2 = 0;
params.k3 = 0;

params.Tfinal = 200;
params.dt = 0.025;
params.Nt = params.Tfinal/params.dt;

params.vP = 0.05;
params.vW = 0.05;


[PSol, WSol] = main(params);
time_vec = (1:params.Nt+1)*params.dt;
Psum = sum(PSol,1)*params.dx;
Wsum = sum(WSol,1)*params.dx;

% damPsum = sum(PSol.*params.x')*params.dx;
% damWsum = sum(WSol.*params.x')*params.dx;
% damPave = damPsum ./ Psum;
% damWave = damWsum ./ Wsum;

f3 = figure(3);
plot(time_vec, log(Psum), 'b-', time_vec, log(Wsum), 'k--');
legend('$\hat{P}$','$\hat{W}$', 'Interpreter', 'latex', 'location', 'best');
xlabel('t', 'Interpreter','latex'); 
ylabel('$\log(\mathrm{Total\,\,number})$', 'Interpreter', 'latex');
title('$\hat{f}=0.20<\tilde{f}$', 'Interpreter','latex');

% pP = polyfit(time_vec(2001:8001), log(Psum(2001:8001)), 1);
% pW = polyfit(time_vec(2001:8001), log(Wsum(2001:8001)), 1);

saveas(f3, 'bifurcation_c', 'fig');
saveas(f3, 'bifurcation_c', 'svg');