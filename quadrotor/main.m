clear all
close all
clc

% adding the subfolders to the path
addpath(genpath('functions'))
addpath(genpath('data'))

% loads:
%    hovering equilibrium (xs,us)
%    continuous time matrices Ac,Bc of the linearization
%    matrices sys.A, sys.B of the inner-loop discretized with sampling period sys.Ts
%    outerController optimizer instance for the outer controller
load('quadData.mat')
outerController = getOuterController(Ac);
disp('Data successfully loaded')

%%%%%%%%%%%%%%%% ADD YOUR CODE BELOW THIS LINE %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define some constants.
[nx, nu] = size(sys.B);                         % State and input dimenstions
T = 10;                                         % Simulation time [s]
bForces = 0;                                    % Determines if FORCES is used

%% %%%%%%%%%%%%%%%%%%%%%    First MPC controller %%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART I - First MPC controller...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization----------------------------------------------------------
x0 = [-1; 0.1745; -0.1745; 0.8727; 0; 0; 0];
N = 10;
Q = eye(nx);
R = eye(nu);
% TUNNING -------
zdot  = 1;
alpha = 9;  % ROLL
beta  = 17; % PTICH
Q(alpha) = 80;
Q(beta) = 10;

% TUNNING -------

% Constraints--------------------------------------------------------------
stateConstraint.Max = [1; 10/180*pi; 10/180*pi; Inf; 15/180*pi; 15/180*pi; 60/180*pi;];
stateConstraint.Min = -[1; 10/180*pi; 10/180*pi; Inf; 15/180*pi; 15/180*pi; 60/180*pi;];
inputConstraint.Max = [1; 1; 1; 1] - us;        % Add Steady State to Constraints
inputConstraint.Min = [0; 0; 0; 0] - us;        % Add Steady State to Constraints
% Create Constraints in Matrix-Form: Hx*x <= kx | Hu*u <= ku
Hx = [eye(nx);-eye(nx)];
kx = [stateConstraint.Max; -stateConstraint.Min];
Hu = [eye(nu);-eye( nu)];
ku = [inputConstraint.Max; -inputConstraint.Min];

% Final Set
[K_inf,P_inf,~] = dlqr(sys.A,sys.B,Q,R);  %Solve Infinite Horizon discrete LQR-Problem
sysCL = LTISystem('A',sys.A - sys.B*K_inf);
P_state = Polyhedron([Hx;-Hu*K_inf],[kx;ku]); % Set stateConstrPolyhedron:  Hx *x <= kx | -Hu*K*x <= ku
sysCL.x.with('setConstraint');
sysCL.x.setConstraint = P_state;
invSet = sysCL.invariantSet();
Hx_f = invSet.A;
kx_f = invSet.b;

% Debugging: Check for Stability:  eig(sys.A - sys.B*K)


%YALMIP
x = sdpvar(nx,N+1);
u = sdpvar(nu,N);
constraints = [];
objective = 0;

for i=1:N
    constraints = [constraints; x(:,i+1) == sys.A*x(:,i) + sys.B*u(:,i)]; % Dynamics
    constraints = [constraints; Hx*x(:,i)<=kx]; % state constraints
    constraints = [constraints; Hu*u(:,i)<=ku]; % input constraints
    
    objective = objective + x(:,i)'*P_inf*x(:,i) + u(:,i)'*R*u(:,i); % costs
end

% Final Constraint / Costs
constraints = [constraints; Hx_f*x(:,N+1)<=kx_f]; % final state constraints
objective = objective + x(:,N+1)'*P_inf*x(:,N+1); % final costs


options = sdpsettings;
innercontroller = optimizer(constraints, objective, options, x(:,1), u(:,1));  % 1st arg: input to funciotn, 2nd arg: output to optimizer;
simQuad(sys, innercontroller, bForces, x0, T);

%% %%%%%%%%%%%%%%%%%%%%%  Reference Tracking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART II - Reference tracking...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Clear variables
clearvars x u objective constraints innerController


% Define parameter for reference tracking
constantReference = 1;
bForces = 1;

% Initialization----------------------------------------------------------
x0 = zeros(nx,1);
N = 30;
Q = eye(nx);
R = eye(nu);
% TUNNING -------
zdot  = 1;
alpha = 9;  % ROLL
beta  = 17; % PTICH
Q(alpha) = 30;
Q(beta)  = 120;

% TUNNING -------

% Final Set
[K_inf,P_inf,~] = dlqr(sys.A,sys.B,Q,R);  %Solve Infinite Horizon discrete LQR-Problem
sysCL = LTISystem('A',sys.A - sys.B*K_inf);
P_state = Polyhedron([Hx;-Hu*K_inf],[kx;ku]); % Set stateConstrPolyhedron:  Hx *x <= kx | -Hu*K*x <= ku
sysCL.x.with('setConstraint');
sysCL.x.setConstraint = P_state;
invSet = sysCL.invariantSet();
Hx_f = invSet.A;
kx_f = invSet.b;


% Constraints--------------------------------------------------------------
stateConstraint.Max = [1; 10/180*pi; 10/180*pi; Inf; 15/180*pi; 15/180*pi; 60/180*pi;];
stateConstraint.Min = -[1; 10/180*pi; 10/180*pi; Inf; 15/180*pi; 15/180*pi; 60/180*pi;];
inputConstraint.Max = [1; 1; 1; 1] - us;        % Add Steady State to Constraints
inputConstraint.Min = [0; 0; 0; 0] - us;        % Add Steady State to Constraints
% Create Constraints in Matrix-Form: Hx*x <= kx | Hu*u <= ku
Hx = [eye(nx);-eye(nx)];
kx = [stateConstraint.Max; -stateConstraint.Min];
Hu = [eye(nu);-eye( nu)];
ku = [inputConstraint.Max; -inputConstraint.Min];


%YALMIP  DELTA FORMULATION

x = sdpvar(nx,N+1);
u = sdpvar(nu,N);
ref = sdpvar(4,1);

% Use formula of L07 slide 14 to compute xr and ur directly, which is not a
% problem, since the matrix has full rank.
C = [eye(4), zeros(4,3)];
Z = [eye(nx) - sys.A, - sys.B; C, zeros(nu)];
ss = Z\[zeros(nx,1); ref];
xr = ss(1:nx);
ur = ss(nx+1:end);      % Nullvektor

constraints = [];
objective = 0;

for i=1:N
    constraints = [constraints; x(:,i+1)-xr == sys.A*(x(:,i)-xr) + sys.B*(u(:,i)-ur)]; % Dynamics
    constraints = [constraints; Hx*x(:,i)<=kx]; % state constraints
    constraints = [constraints; Hu*u(:,i)<=ku]; % input constraints
    
    objective = objective + (x(:,i)-xr)'*P_inf*(x(:,i)-xr) + (u(:,i)-ur)'*R*(u(:,i)-ur); % costs
end

% Final Constraint / Costs
constraints = [constraints; Hx_f*x(:,N+1)<=kx_f]; % final state constraints
objective = objective + (x(:,N+1)-xr)'*P_inf*(x(:,N+1)-xr); % final costs


% Simulate either with constant or varying reference
x0 = zeros(nx,1);
if constantReference
    r = [1.0; 0.1745; -0.1745; 1.7453];
else
    T_vec = 0:sys.Ts:T;
    r = [1.0*ones(size(T_vec)); 0.1745*sin(T_vec); -0.1745*sin(T_vec); pi/2*ones(size(T_vec))];
end

%% 
clearvars  innerController
options = sdpsettings;
innerController = optimizer(constraints, objective, options, [x(:,1)', ref']', u(:,1))  % 1st arg: input to funciotn, 2nd arg: output to optimizer;


% Create controller object (generates code)
% for a complete list of codeoptions, see 
% https://www.embotech.com/FORCES-Pro/User-Manual/Low-level-Interface/Solver-Options
codeoptions = getOptions('simpleMPC_solver'); % give solver a name
% innerController = optimizerFORCES(constraints, objective, codeoptions, [x(:,1)', ref']', u(:,1), {'xinit'}, {'u0'});
% [output,exitflag,info] = simpleMPC_solver({[x0',r']'});
tic
simpleMPC_solver({[x0',r']'})
toc
%vs
tic
innerController({[x0',r']'})
toc

tic
simQuad(sys, innerController, bForces, x0, T, r);
toc
%% %%%%%%%%%%%%%%%  First simulation of the nonlinear model %%%%%%%%%%%%%%%%%
fprintf('PART III - First simulation of the nonlinear model...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run the nonlinear simulation
sim('simulation1')

%% %%%%%%%%%%%%%%%%%%%%%%%  Offset free MPC  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART IV - Offset free MPC...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%  Simulation of the nonlinear model %%%%%%%%%%%%%%%%%%%%
fprintf('PART V - simulation of the nonlinear model...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%  Slew Rate Constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART VI - Slew Rate Constraints...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%  Soft Constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART VII - Soft Constraints...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FORCES Pro %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART VIII - FORCES Pro...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xt ut t rt deltat] 

