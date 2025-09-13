%% run_calibrated_sims.m

clear; close all; clc;

params = initializeParams;

% set simulation choices to the values specified
params.Nx = 400;
params.x = linspace(0,2,params.Nx+1);
params.dx = params.x(2) - params.x(1);
params.Tfinal = 400;
params.dt = 0.025;
params.Nt = round(params.Tfinal / params.dt);
params.vP = 0.05;
params.vW = 0.05;
params.lambdaR = 0;
params.delta = 0.6 * params.x';
params.k3 = 0;
params.k4 = 0;
params.k1 = 1.4087e-3;
params.k2 = 0.1 * params.k1;
params.lambdaP = 1;
params.p1 = 0.5;
params.p2 = 1 - params.p1;
params.p3 = 0;

% targets
targetW = 10;
targetRatio = 7;

case_list= 1:5;

results = struct();

for idx = case_list
    fprintf('------Running BC case %d -----\n', idx);
    locParams = params;

    [LocParams, calib_info] = calibrate_k1_lambdaP(locParams, idx, targetW, targetRatio);
    fprintf('Calibration finished: k1 = %.4f, lambdaP= %.4f (iters=%d)\n', locParams.k1, locParams.lambdaP, calib_info.iter);

    [PSol, WSol] = main(locParams, idx);
    Psum = sum(PSol(:,end)) * locParams.dx;
    Wsum = sum(WSol(:,end)) * locParams.dx;
    ratio = Wsum / Psum;

    results(idx).params = PSol;
    results(idx).PSol = PSol;
    results(idx).WSol = WSol;
    results(idx).Psum = Psum;
    results(idx).Wsum = Wsum;
    results(idx).ratio = ratio;
    results(idx).calib = calib_info;

    fprintf('After full run: Psum = %.4f, Wsum = %.4f, ratio = %.4f', Psum, Wsum, ratio);
end

save('calibrated_runs.mat', 'results');
