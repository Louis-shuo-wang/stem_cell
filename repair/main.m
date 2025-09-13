function [PSol, WSol] = main(params)

if nargin == 0
    params = initializeParams;
end

% initial condition
x = params.x;
Nx = params.Nx;

P0 = zeros(size(x));
W0 = zeros(size(x));
P0(floor(Nx/10)+1:2*floor(Nx/10))=10;
W0(floor(Nx/10)+1:2*floor(Nx/10))=10;


PSol = zeros(params.Nx+1, params.Nt+1);
WSol = zeros(params.Nx+1, params.Nt+1);
PSol(:,1) = P0;
WSol(:,1) = W0;


for t = 1:params.Nt
    Y = [PSol(:,t);WSol(:,t)];
    dYdt = rhs_massconserve(t,Y,params);
    PSol(:,t+1) = PSol(:,t) + params.dt * dYdt(1:Nx+1);
    WSol(:,t+1) = WSol(:,t) + params.dt * dYdt(Nx+2:end);
    PSol(:,t+1) = max(0, PSol(:,t+1));
    WSol(:,t+1) = max(0, WSol(:,t+1));
end

end

function dYdt = rhs_massconserve(t,Y,params)
Nx = params.Nx;
dx = params.dx;
x = params.x;
vP = params.vP;
vW = params.vW;
lambdaP = params.lambdaP;
k1 = params.k1;
k2 = params.k2;
k3 = params.k3;
k4 = params.k4;
m1 = params.m1;
m2 = params.m2;
m3 = params.m3;
m4 = params.m4;

p1=params.p1;
p2=params.p2;
p3=params.p3;
alpha1=params.alpha1;
alpha2=params.alpha2;
beta1=params.beta1;
beta2=params.beta2;
gamma1=params.gamma1;
gamma2=params.gamma2;
lambdaR=params.lambdaR;
delta=params.delta;
rho=params.rho;


% unpack state vector
P = Y(1:Nx+1);
W = Y(Nx+2:end);
Psum = sum(P) * dx;
Wsum = sum(W) * dx;
p1 = p1/(1+(k1*Wsum)^m1);
p2 = p2/(1+(k2*Wsum)^m2);
p3 = 1-p1-p2;
lambdaP = lambdaP/(1+(k3*Wsum)^m3);
lambdaR = lambdaR/(1+(k4*Psum)^m4);
%--------------------------------------
% Conservative advection (upwind flux)
%--------------------------------------
% --- P transport ---
fluxP = zeros(Nx+2,1);  % define at interfaces (Nx+2 points)
if vP >= 0
    fluxP(2:Nx+1) = vP * P(1:Nx);    % upwind: left state
else
    fluxP(2:Nx+1) = vP * P(2:Nx+1);  % upwind: right state
end
% divergence of flux (Nx+1 points)
divFluxP = (fluxP(2:Nx+2) - fluxP(1:Nx+1)) / dx;

rhsP = -lambdaP*P - divFluxP;

% --- W transport ---
fluxW = zeros(Nx+2,1);
if vW >= 0
    fluxW(2:Nx+1) = vW * W(1:Nx);
else
    fluxW(2:Nx+1) = vW * W(2:Nx+1);
end
divFluxW = (fluxW(2:Nx+2) - fluxW(1:Nx+1)) / dx;

rhsW = - delta.*W - divFluxW;

%------------------------------------------------
% Nonlocal degradation repair during rejuvenation
%------------------------------------------------

rhsP = rhsP + lambdaR * scaling(W, x, rho);
rhsW = rhsW - lambdaR * scaling(W, x, rho);

%--------------------------------------
% Nonlocal proliferation / maturation
%--------------------------------------
% Symmetric self-renewal (p1)
rhsP = rhsP + p1*lambdaP * scaling(P, x, alpha1);
rhsP = rhsP + p1*lambdaP * scaling(P, x, alpha2);

% Symmetric maturation (p2 -> W)
rhsW = rhsW + p2*lambdaP * scaling(P, x, beta1);
rhsW = rhsW + p2*lambdaP * scaling(P, x, beta2);

% Asymmetric division (p3)
rhsP = rhsP + p3*lambdaP * scaling(P, x, gamma1);
rhsW = rhsW + p3*lambdaP * scaling(P, x, gamma2);

% pack
dYdt = [rhsP; rhsW];
end

function term = scaling(P, x, alpha)
Nx = length(x);
dx = x(2)-x(1);
term = zeros(size(P));
for j = 1:Nx
    x_new = alpha * x(j);
    if x_new >= x(end), continue; end
    idx = floor(x_new/dx) + 1; % left bin (1-based index in Matlab)
    idx = min(idx, Nx); % avoid overflow when idx+1 > Nx+1
    w = (x_new - (idx-1)*dx)/dx; % interpolation weight
    term(idx) = term(idx) + (1-w)*P(j);
    term(idx+1) = term(idx+1) + w*P(j);
end
end