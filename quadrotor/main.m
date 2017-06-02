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

%%  ADD YOUR CODE BELOW THIS LINE %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define some constants.
[nx, nu] = size(sys.B); % State and input dimenstions
T = 10;                 % Simulation time [s]
N = 20;                 % Time horizon

% Define cost parameters
Q = diag([1000 8000 6000 30 4 10 10]);
R = 0.1*eye(nu);

% R(2,2) = 0.01;
% R(3,3)  = 0.01;

% Define constraints
stateConstraint.Max = [1; 10/180*pi; 10/180*pi; Inf; 15/180*pi; 15/180*pi; 60/180*pi;];
stateConstraint.Min = -[1; 10/180*pi; 10/180*pi; Inf; 15/180*pi; 15/180*pi; 60/180*pi;];
inputConstraint.Max = [1; 1; 1; 1] - us;        % Add Steady State to Constraints
inputConstraint.Min = [0; 0; 0; 0] - us;        % Add Steady State to Constraints

% Create Constraints in Matrix-Form: Hx*x <= kx | Hu*u <= ku
Hx = [eye(nx);-eye(nx)];
kx = [stateConstraint.Max; -stateConstraint.Min];
Hu = [eye(nu);-eye( nu)];
ku = [inputConstraint.Max; -inputConstraint.Min];

% Compute the final set
[K_inf,P_inf,~] = dlqr(sys.A,sys.B,Q,R);  % Solve Infinite Horizon discrete LQR-Problem
sysCL = LTISystem('A',sys.A - sys.B*K_inf);
P_state = Polyhedron([Hx;-Hu*K_inf],[kx;ku]); % Set stateConstrPolyhedron:  Hx *x <= kx | -Hu*K*x <= ku
sysCL.x.with('setConstraint');
sysCL.x.setConstraint = P_state;
invSet = sysCL.invariantSet();
Hx_f = invSet.A;
kx_f = invSet.b;

 %%%%%%%%%%%%%%%%%%%%%    First MPC controller %%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART I - First MPC controller...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define yalmip variables for states and inputs
x = sdpvar(nx,N+1);
u = sdpvar(nu,N);

% Init constraints and objective function
constraints = [];objective = 0;
for i=1:N
    % Add system dynamics
    constraints = [constraints; x(:,i+1) == sys.A*x(:,i) + sys.B*u(:,i)];
    % Add state constraints
    constraints = [constraints; Hx*x(:,i)<=kx];
    % Add input constraints
    constraints = [constraints; Hu*u(:,i)<=ku];
    % Add cost
    objective = objective + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i);
end
% Add final constraints
constraints = [constraints; Hx_f*x(:,N+1)<=kx_f];
% Add final cost
objective = objective + x(:,N+1)'*P_inf*x(:,N+1);

% Construct optimizer
options = sdpsettings;
innercontroller = optimizer(constraints, objective, options, x(:,1), u(:,1));

% Simulate with initial condition
x0 = [-1; 0.1745; -0.1745; 0.8727; 0; 0; 0];
simQuad(sys, innercontroller, 0, x0, T);

%% %%%%%%%%%%%%%%%%%%%%%  Reference Tracking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART II - Reference tracking...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear variables
clearvars x u objective constraints innerController



% Define yalmip states, inputs and reference
x = sdpvar(nx,N+1);
u = sdpvar(nu,N);
ref = sdpvar(4,1);

% Use formula of L07 slide 14 to compute xr and ur directly, which is not a
% problem, since the matrix has full rank.
C = [eye(4), zeros(4,3)];
Z = [eye(nx) - sys.A, - sys.B; C, zeros(nu)];
ss = Z\[zeros(nx,1); ref];
xr = ss(1:nx);
ur = ss(nx+1:end);

% Init constraints and objective function
constraints = [];
objective = 0;
for i=1:N
    % Add system dynamics
    constraints = [constraints; x(:,i+1)-xr == sys.A*(x(:,i)-xr) + sys.B*(u(:,i)-ur)];
    % Add state constraints
    constraints = [constraints; Hx*x(:,i)<=kx];
    % Add input constraints
    constraints = [constraints; Hu*u(:,i)<=ku];
    % Add state cost
    objective = objective + (x(:,i)-xr)'*Q*(x(:,i)-xr) + (u(:,i)-ur)'*R*(u(:,i)-ur);
end
% Add final constraints and cost
constraints = [constraints; Hx*x(:,N+1)<=kx];
objective = objective + (x(:,N+1)-xr)'*P_inf*(x(:,N+1)-xr);

% Construct optimizer
options = sdpsettings;
innerController = optimizer(constraints, objective, options, [x(:,1)', ref']', u(:,1));

% Simulate with constant and varying reference
x0 = zeros(nx,1);
r_const = [1.0; 0.1745; -0.1745; 1.7453];
simQuad(sys, innerController, 0, x0, T, r_const);
T_vec = 0:sys.Ts:T;
r_var = [1.0*ones(size(T_vec)); 0.1745*sin(T_vec); -0.1745*sin(T_vec); pi/2*ones(size(T_vec))];
simQuad(sys, innerController, 0, x0, T, r_var);

%% %%%%%%%%%%%%%%%  First simulation of the nonlinear model %%%%%%%%%%%%%%%
fprintf('PART III - First simulation of the nonlinear model...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulate the nonlinear model
sim('simulation1')

%% %%%%%%%%%%%%%%%%%%%%%  Offset free MPC  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART IV - Offset free MPC...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars x u objective constraints innerController ref ss xr ur

seed = 2;
rng(seed);

% Define variables for controller
ref = sdpvar(4,1);      % reference
x = sdpvar(nx,N+1);     % states
u = sdpvar(nu,N);       % inputs
d_est = sdpvar(nx,1);       % disturbance

% Define system equations as x(k+1) = A x(k) + B u(k) + B_d d(k), 
% y(k) = C x(k) + C_d d(k) and the disturbance dynamics as d(k+1) = d(k).
C = eye(nx);
B_d = eye(nx);
C_d = zeros(nx);

% Define state observer dynamics as
% [x_hat(k+1); d_hat(k+1)] = (A_aug - L C_aug) [x_hat(k); d_hat(k)] +
% [B_aug L] [u(k) x(k)], where L is chosen such thath the state observer
% dynamics are stable and go to zero.
A_aug = [sys.A B_d; zeros(nx) eye(nx)];
B_aug = [sys.B; zeros(nx,nu)];
C_aug = [eye(nx) eye(nx)];

Q_ = diag([1*ones(1,nx) [50 1 1 500 10 10 0.01 ]]);
% Q(1,1) = 10;
% Q(3,3) = 100;
% Q(4,4) = 5;
R_ = 10*eye(nx);
% R(1,1) = 0.10;
% R(2,2) = 0.001;
% R(3,3)  = 0.001;
% R(4,4) = 0.01;

L = dlqr(A_aug',C_aug',Q_,R_)';
% L = [eye(nx); diag([0.1 0.1 1 0.1 1 1 1])];

% p = [0.91, 0.92, 0.93, 0.901, 0.88, 0.8, 0.99, 0.95, 0.955, 0.99, 0.99, 0.98, 0.97, 0.99]';
% L = place(A_aug',C_aug',p)'

% Defining the filter
filter.Af = A_aug - L*C_aug;
abs(eig(filter.Af))
filter.Bf = [B_aug L];

% Use formula of L07 slide 29 to compute xr and ur directly, which is not a
% problem, since the matrix has full rank.
Z = [sys.A-eye(nx), sys.B; C(1:4,:), zeros(nu)];
ss = Z\[-B_d*d_est; ref - C_d(1:4,:)*d_est];
xr = ss(1:nx);
ur = ss(nx+1:end);

% Init constraints and objective function
constraints = [];
objective = 0;
for i = 1:N
    % Add state evolution constraints.
    constraints = [constraints, x(:,i+1) == sys.A*(x(:,i)) + sys.B*(u(:,i)) + B_d*d_est];
    % Add state constraints.
    constraints = [constraints, Hx*(x(:,i)) <= kx];
    % Add input constraints.
    constraints = [constraints, Hu*(u(:,i)) <= ku];
    
    % Add to objective function
    objective = objective + (x(:,i)-xr)' * Q * (x(:,i)-xr) + (u(:,i)-ur)' * R * (u(:,i)-ur);
end
% Add final constraints and objective
constraints = [constraints, Hx*(x(:,N+1)) <= kx];
objective = objective + (x(:,i)-xr)' * P_inf * (x(:,i)-xr);

% Construct optimizer
options = sdpsettings;
innerController = optimizer(constraints, objective, options, [x(:,1)', ref', d_est']', u(:,1));

% Simulate either with constant or varying reference
x0 = zeros(nx,1);
r_const = [0.8; 0.12; -0.12; pi/2];
simQuad(sys, innerController, 0, x0, T, r_const, filter, []);
T_vec = 0:sys.Ts:T;
r_var = [0.8*ones(size(T_vec)); 0.12*sin(T_vec); -0.12*sin(T_vec); pi/2*ones(size(T_vec))];
simQuad(sys, innerController, 0, x0, T, r_var, filter, []);

%% %%%%%%%%%%%%%%%%  Simulation of the nonlinear model %%%%%%%%%%%%%%%%%%%%
fprintf('PART V - simulation of the nonlinear model...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulate the nonlinear model
sim('simulation2')

%% %%%%%%%%%%%%%%%%%%%%%  Slew Rate Constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART VI - Slew Rate Constraints...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%  Soft Constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART VII - Soft Constraints...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%  FORCES Pro %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART VIII - FORCES Pro...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create controller object (generates code)
% for a complete list of codeoptions, see 
% https://www.embotech.com/FORCES-Pro/User-Manual/Low-level-Interface/Solver-Options
% codeoptions = getOptions('simpleMPC_solver'); % give solver a name
% innerController = optimizerFORCES(constraints, objective, codeoptions, [x(:,1)', ref']', u(:,1), {'xinit'}, {'u0'});
% [output,exitflag,info] = simpleMPC_solver({[x0',r']'});
% [xt ut t rt deltat] x
