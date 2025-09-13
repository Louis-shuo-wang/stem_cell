function params = initializeParams()

% spatial discretization
params.Nx = 200;
params.x = linspace(0,2,params.Nx+1);
params.dx = params.x(2)-params.x(1);

% temporal discretization
params.Tfinal = 500;
params.dt = 0.1;
params.Nt=params.Tfinal/params.dt;
params.CFL = 0.5;

params.vP = 0.025; 
params.vW = 0.025;
params.dt = min(params.dt, params.CFL*params.dx/max(params.vP,params.vW));
params.k1 = 0;
params.k2 = 0;
params.k3 = 0;
params.k4 = 0;
params.m1 = 2;
params.m2 = 2;
params.m3 = 2;
params.m4 = 2;

params.alpha1=0.5; 
params.alpha2=0.5;
params.beta1=0.5; 
params.beta2=0.5;
params.gamma1=1/3; 
params.gamma2=2/3;

params.lambdaP = 1;
params.lambdaR = 0;          % could be function of x
params.delta   = 0;              % increasing with x


params.p1=0.5; 
params.p2=0.5; 
params.p3=0;  % base case