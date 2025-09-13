% This matlab function is used to plot degradation-inducing effect of
% rejuvenation

clear all; close all;
params = initializeParams;
params.Nx = 400;
params.x = linspace(0,2,params.Nx+1);
params.dx = params.x(2) - params.x(1);

params.Tfinal = 800;
params.dt = 0.025;
params.Nt = params.Tfinal/params.dt;

params.vP = 0.05;
params.vW = 0.05;

params.lambdaR = 0;
params.k3 = 0;
params.k4 = 0;

i = 2;


params.p1 = 0.3;
params.p2 = 1-params.p1;
params.p3 = 0;
params.lambdaP = 1;
params.lambdaR = 0;
params.delta = (1 + 10*params.lambdaP/7)/10;
params.k1 = 1/(10*sqrt(7));
params.k2 = params.k1;

[PSol,WSol] = main(params, i);
Psum = sum(PSol,1)*params.dx;
Wsum = sum(WSol,1)*params.dx;
ratio = Wsum./Psum;

fprintf("ratio is %.4f\n", ratio(end));
fprintf('W number is %.4f\n', Wsum(end));
fprintf('P number is %.4f\n', Psum(end));
fprintf('k1 should be %.4e\n', Wsum(end)*params.k1/10);

% lambdaPs = [1.084875, 1.08775, 1.09016, 1.09243, 1.09459];
% lambdaRs = [0.01, 0.03, 0.05, 0.07, 0.09];
% p1s = [0.5, 0.5, 0.5, 0.5, 0.5];
% k1s = [0.00003888, 0.00006577, 0.000084438, 0.000099331, 0.0001119615];


for t = 1:params.Nt/100:params.Nt+1
    plot(params.x,PSol(:, t));
    pause(0.05);
end
