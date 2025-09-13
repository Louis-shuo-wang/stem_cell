clear all;
close all;
% parameters
hatlambdaP = 0.4;
lambdaR = 0.02;
delta = 0.2;
hatp = 0.25;

k1 = 2;
k2 = 2;
m1 = 2;
m2 = 2;

Tfinal = 50;
dt = 0.05;
N = Tfinal/dt;

P0 = 2.5;
W0 = 1;
P = zeros(N+1,1);
W = zeros(N+1,1);

P(1) = P0;
W(1) = W0;

for i = 1:N
p = (2*hatp-1)/(1+(k1*W(i))^m1);
lambdaP = hatlambdaP/(1+(k2*W(i))^m2);
    P(i+1) = P(i) + dt*( p*lambdaP * P(i) + lambdaR * W(i));
    W(i+1) = W(i) + dt*( (1-p)*lambdaP*P(i) - lambdaR*W(i) - delta * W(i));
end

time_vec = (0:N)*dt;
plot(time_vec, P, 'b-', time_vec, W, 'k--');
xlabel('t', 'FontSize', 15);
ylabel('Number', 'FontSize', 15);
ylim([0 6]);
legend('$\hat{P}$', '$\hat{W}$', 'FontSize', 15);
title('Total number', 'FontSize', 16);