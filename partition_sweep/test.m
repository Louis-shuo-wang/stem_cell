% This matlab function is used to plot degradation-inducing effect of
% rejuvenation

clear all; close all;
params = initializeParams;
params.Nx = 400;
params.x = linspace(0,2,params.Nx+1);
params.dx = params.x(2) - params.x(1);

params.Tfinal = 2000;
params.dt = 0.025;
params.Nt = params.Tfinal/params.dt;

params.vP = 0.02;
params.vW = 0.02;

params.delta = 0.6 * params.x';
params.k3 = 0;
params.k4 = 0;

parama.p1 = 0.5;
params.p2 = 0.5;
params.p3 = 0;

params.lambdaP = 0.6825;

params.lambdaR = 0;

params.alpha1=0.5; 
params.beta1=0.3; 
params.gamma1=1/3; 
params.alpha2=1-params.alpha1;
params.beta2=1-params.beta1;
params.gamma2=1-params.gamma1;

params.k1 = 9.5442e-06;
params.k2 = 0.1*params.k1;

[PSol,WSol] = main(params);
Psum = sum(PSol,1)*params.dx;
Wsum = sum(WSol,1)*params.dx;
ratio = Wsum./Psum;

fprintf("ratio is %.4f\n", ratio(end));
fprintf('W number is %.4f\n', Wsum(end));
fprintf('P number is %.4f\n', Psum(end));
fprintf('k1 should be %.4e\n', Wsum(end)*params.k1/10);

% for t = 1:10:params.Nt+1
%     plot(params.x,WSol(:,t));
%     pause(0.05);
% end