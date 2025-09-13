%% the right panel
clear; close all; 
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

params.Tfinal = 100;
params.dt = 0.025;
params.Nt = params.Tfinal/params.dt;

params.vP = 0.2;
params.vW = 0.2;

P1Sol = zeros(params.Nx+1,params.Nt+1);
W1Sol = zeros(params.Nx+1,params.Nt+1);

P0 = zeros(params.Nx+1,1);
W0 = zeros(params.Nx+1,1);
P0(floor(params.Nx/10)+1:2*floor(params.Nx/10))=10;
W0(floor(params.Nx/10)+1:2*floor(params.Nx/10))=10;

P1Sol(:,1) = P0;
W1Sol(:,1) = W0;

for t = 1:params.Nt
P = P1Sol(:, t);
W = W1Sol(:, t);

fluxP = zeros(params.Nx+2,1); % define at interfaces (Nx+2 points)
if params.vP >= 0
fluxP(2:params.Nx+1) = params.vP * P(1:params.Nx); % upwind: left state
else
fluxP(2:params.Nx+1) = params.vP * P(2:params.Nx+1); % upwind: right state
end
% divergence of flux (Nx+1 points)
divFluxP = (fluxP(2:params.Nx+2) - fluxP(1:params.Nx+1)) / params.dx;

% --- W transport ---
fluxW = zeros(params.Nx+2,1);
if params.vW >= 0
fluxW(2:params.Nx+1) = params.vW * W(1:params.Nx);
else
fluxW(2:params.Nx+1) = params.vW * W(2:params.Nx+1);
end
divFluxW = (fluxW(2:params.Nx+2) - fluxW(1:params.Nx+1)) / params.dx;

rhsP = - divFluxP + (2*params.p1-1)*params.lambdaP*P1Sol(:,t)+params.lambdaR*W1Sol(:,t);
rhsW = - divFluxW + (2-2*params.p1)*params.lambdaP*P1Sol(:,t)-(params.delta+params.lambdaR)*W1Sol(:,t);

P1Sol(:,t+1) = P1Sol(:,t) + params.dt*rhsP;
W1Sol(:,t+1) = W1Sol(:,t) + params.dt*rhsW;
end

P1sum = sum(P1Sol,1)*params.dx;
W1sum = sum(W1Sol,1)*params.dx;

time_vec = (1:params.Nt+1)*params.dt;

damP1sum = sum(P1Sol.*params.x')*params.dx;
damW1sum = sum(W1Sol.*params.x')*params.dx;
damP1ave = damP1sum ./ P1sum;
damW1ave = damW1sum ./ W1sum;

set(groot, 'defaultAxesFontSize', 30);
set(groot, 'defaultTextFontSize', 30);
set(groot, 'defaultAxesLineWidth', 1.2);
set(groot, 'defaultLineLinewidth', 2.5);
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesBox', 'off');

% plot(time_vec, damP1ave, 'b-', time_vec, damW1ave, 'k--');
% legend('dam_P_ave','dam_W_ave');
% 
% pP = polyfit(time_vec(1:60/params.dt), damP1ave(1:60/params.dt), 1);   
% pW = polyfit(time_vec(1:60/params.dt), damW1ave(1:60/params.dt), 1); 

legendLabels = {};
times = [10,35,60]/params.dt;
colors = lines(length(times));

figure;

for i = 1:length(times)
    time = times(i);
    plot(params.x,P1Sol(:,time), '-','Color',colors(i,:),'LineWidth',2); hold on;
    legendLabels{end+1} = sprintf('P at t= %.0f', time*params.dt); 

    plot(params.x,W1Sol(:,time), '--','Color',colors(i,:),'LineWidth',2);
    legendLabels{end+1} = sprintf('W at t= %.0f', time*params.dt);
end

    xlabel('x','Interpreter','latex');
    ylabel('Density', 'Interpreter','latex');
    title('Damage Distribution', 'Interpreter','latex');
    legend(legendLabels, 'Interpreter','latex', 'location','best');

saveas(gcf, 'bounded_b.fig');
saveas(gcf, 'bounded_b.svg');

% After loop: show space-time images
% figure;
% for t = 1:20:params.Nt+1
% subplot(1,2,1);
% plot(params.x, P1Sol(:,t));
% axis xy; xlabel('x'); ylabel('P'); title(sprintf('P(x,%.0f)', (t-1)*params.dt));
% 
% subplot(1,2,2);
% plot(params.x, W1Sol(:,t));
% axis xy; xlabel('x'); ylabel('W'); title(sprintf('W(x,%.0f)', (t-1)*params.dt));
% pause(0.05);
% end
