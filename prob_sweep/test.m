% This matlab function is used to plot degradation-inducing effect of
% rejuvenation

clear all; close all;
params = initializeParams;
params.Nx = 400;
params.x = linspace(0,2,params.Nx+1);
params.dx = params.x(2) - params.x(1);

params.Tfinal = 400;
params.dt = 0.025;
params.Nt = params.Tfinal/params.dt;

params.vP = 0.05;
params.vW = 0.05;

params.delta = 0.6 * params.x';
params.k3 = 0;
params.k4 = 0;

params.lambdaP = 1.08357;
% in this threshold sweep, we change the threshold value and fix
% \lambda_R=0.05
params.lambdaR = 0;


params.p3 = 0.9;
parama.p1 = (1-params.p3)/2;
params.p2 = (1-params.p3)/2;
params.k1 =  0.00158114;
params.k2 = 0.1*params.k1;

[PSol,WSol] = main(params);
Psum = sum(PSol,1)*params.dx;
Wsum = sum(WSol,1)*params.dx;
ratio = Wsum./Psum;

fprintf("ratio is %.4f\n", ratio(end));
fprintf('W number is %.4f\n', Wsum(end));
fprintf('P number is %.4f\n', Psum(end));

% lambdaPs = [1.084875, 1.08775, 1.09016, 1.09243, 1.09459];
% lambdaRs = [0.01, 0.03, 0.05, 0.07, 0.09];
% p1s = [0.5, 0.5, 0.5, 0.5, 0.5];
% k1s = [0.00003888, 0.00006577, 0.000084438, 0.000099331, 0.0001119615];

