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

%% %%%%%%%%%%%%%%% ADD YOUR CODE BELOW THIS LINE %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define some constants.
[nx, nu] = size(sys.B); % State and input dimenstions
T = 10;                 % Simulation time [s]
N = 30;                 % Time horizon

% Define cost parameters
Q = diag([10000 4000 9000 5 250 200 1]);
R = 0.1*eye(nu);

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

%% %%%%%%%%%%%%%%%%%%%%%    First MPC controller %%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART I - First MPC controller...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define yalmip variables for x and u
x = sdpvar(nx,N+1);
u = sdpvar(nu,N);

% Init constraints and objective function
constraints = [];
objective = 0;
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
innerController = optimizer(constraints, objective, options, x(:,1), u(:,1));

% Simulate with initial condition
x0 = [-1; 0.1745; -0.1745; 0.8727; 0; 0; 0];
simQuad(sys, innerController, 0, x0, T);

%% %%%%%%%%%%%%%%%%%%%%%  Reference Tracking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART II - Reference tracking...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear variables
clearvars x u objective constraints innerController

% Define yalmip x, u and reference
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
[~, ~, ~, rt_task5 ,~] = simQuad(sys, innerController, 0, x0, T, r_const)
T_vec = 0:sys.Ts:T;
r_var = [1.0*ones(size(T_vec)); 0.1745*sin(T_vec); -0.1745*sin(T_vec); pi/2*ones(size(T_vec))];
simQuad(sys, innerController, 0, x0, T, r_var);

% Construct FORCES Pro optimizer
codeoptions = getOptions('simpleMPC_solver'); % give solver a name
innerController_FORCES = optimizerFORCES(constraints, objective, codeoptions, [x(:,1)', ref']', u(:,1));%, {'xinit'}, {'u0'});

%% %%%%%%%%%%%%%%%  First simulation of the nonlinear model %%%%%%%%%%%%%%%
fprintf('PART III - First simulation of the nonlinear model...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulate the nonlinear model
sim('simulation1')

%% %%%%%%%%%%%%%%%%%%%%%  Offset free MPC  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART IV - Offset free MPC...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear previous variables.
clearvars x u objective constraints innerController ref ss xr ur

% Define variables for controller
ref = sdpvar(4,1);      % reference
x = sdpvar(nx,N+1);     % x
u = sdpvar(nu,N);       % u
dest = sdpvar(nx,1);   % disturbance

% Define system equations as x(k+1) = A x(k) + B u(k) + B_d d(k), 
% y(k) = C x(k) + C_d d(k) and the disturbance dynamics as d(k+1) = d(k).
C = eye(nx);
B_d = eye(nx);
C_d = zeros(nx);        % Set to Zero!!

% Define state observer dynamics as
% [x_hat(k+1); d_hat(k+1)] = (A_aug - L C_aug) [x_hat(k); d_hat(k)] +
% [B_aug L] [u(k) x(k)], where L is chosen such thath the state observer
% dynamics are stable and go to zero.
A_aug = [sys.A B_d; zeros(nx) eye(nx)];
B_aug = [sys.B; zeros(nx,nu)];
C_aug = [eye(nx) eye(nx)];

% Choose cost matrices for observer design
Q_ = diag([1*ones(1,nx) [50 1 1 500 10 10 0.01 ]]);
R_ = 10*eye(nx);

% For this Q_ and R_ also the offset-free trackign with r_const is feasible
% Q_ = diag(ones(2*nx,1));
% R_ = diag(ones(nx,1));
L = dlqr(A_aug',C_aug',Q_,R_)';

% Defining the filter
filter.Af = A_aug - L*C_aug;
% abs(eig(filter.Af))
filter.Bf = [B_aug L];

% Use formula of L07 slide 29 to compute xr and ur directly, which is not a
% problem, since the matrix has full rank.
Z = [sys.A-eye(nx), sys.B; C(1:4,:), zeros(nu)];
ss = Z\[-B_d*dest; ref - 0*C_d(1:4,:)*dest];
xr = ss(1:nx);
ur = ss(nx+1:end);

% Init constraints and objective function
constraints = [];
objective = 0;
for i = 1:N
    % Add state evolution constraints.
    constraints = [constraints, x(:,i+1) == sys.A*(x(:,i)) + sys.B*(u(:,i)) + B_d*dest];
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
innerController = optimizer(constraints, objective, options, [x(:,1)', ref', dest']', u(:,1));

% Simulate either with constant or varying reference
x0 = zeros(nx,1);
r_const = [0.8; 0.12; -0.12; pi/2];
[~, ~, ~, rt_task9 ,~] = simQuad(sys, innerController, 0, x0, T, r_const, filter, [])
T_vec = 0:sys.Ts:T;
r_var = [0.8*ones(size(T_vec)); 0.12*sin(T_vec); -0.12*sin(T_vec); pi/2*ones(size(T_vec))];
simQuad(sys, innerController, 0, x0, T, r_var, filter, []);

% Construct FORCES Pro optimizer
codeoptions = getOptions('offsetFreeMPC_solver'); % give solver a name
innerController_FORCES_task9 = optimizerFORCES(constraints, objective, codeoptions, [x(:,1)', ref', dest']', u(:,1));%, {'xinit'}, {'u0'});

%% %%%%%%%%%%%%%%%%  Simulation of the nonlinear model %%%%%%%%%%%%%%%%%%%%
fprintf('PART V - simulation of the nonlinear model...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulate the nonlinear model
sim('simulation2')

%% %%%%%%%%%%%%%%%%%%%%%  Slew Rate Constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART VI - Slew Rate Constraints...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear previous variables.
clearvars x u objective constraints innerController ref ss xr ur dest

% Define variables for controller
ref = sdpvar(4,1);      % reference
x = sdpvar(nx,N+1);     % x
u = sdpvar(nu,N);       % u
dest = sdpvar(nx,1);   % disturbance
u_prev = sdpvar(nu,1);  % previous input

% Define system equations as x(k+1) = A x(k) + B u(k) + B_d d(k), 
% y(k) = C x(k) + C_d d(k) and the disturbance dynamics as d(k+1) = d(k).
C = eye(nx);
B_d = eye(nx);
C_d = zeros(nx);        % Set to Zero!!

% Define state observer dynamics as
% [x_hat(k+1); d_hat(k+1)] = (A_aug - L C_aug) [x_hat(k); d_hat(k)] +
% [B_aug L] [u(k) x(k)], where L is chosen such thath the state observer
% dynamics are stable and go to zero.
A_aug = [sys.A B_d; zeros(nx) eye(nx)];
B_aug = [sys.B; zeros(nx,nu)];
C_aug = [eye(nx) eye(nx)];

% Choose cost matrices for observer design
Q_ = diag([1*ones(1,nx) [50 1 1 500 10 10 0.01 ]]);
R_ = 10*eye(nx);

% For this Q_ and R_ also the offset-free trackign with r_const is feasible
% Q_ = diag(ones(2*nx,1));
% R_ = diag(ones(nx,1));
L = dlqr(A_aug',C_aug',Q_,R_)';

% Defining the filter
filter.Af = A_aug - L*C_aug;
% abs(eig(filter.Af))
filter.Bf = [B_aug L];

% Use formula of L07 slide 29 to compute xr and ur directly, which is not a
% problem, since the matrix has full rank.
Z = [sys.A-eye(nx), sys.B; C(1:4,:), zeros(nu)];
ss = Z\[-B_d*dest; ref - 0*C_d(1:4,:)*dest];
xr = ss(1:nx);
ur = ss(nx+1:end);

% Init constraints and objective function
constraints = [];
objective = 0;
for i = 1:N
    % Add state evolution constraints.
    constraints = [constraints, x(:,i+1) == sys.A*(x(:,i)) + sys.B*(u(:,i)) + B_d*dest];
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

% Add slew rate constraints to previously defined contraints.  For delta =
% 0.25125, roll and pitch angles start to overshoot and for delta = 0.25
% the problem is infeasible. delta = 0.26 is okay.
delta = 0.26*ones(nu,1);
constraints = [constraints, (-delta <= u_prev - u(:,1) <= delta):'slew rate'];
for i = 2:N
    constraints = [constraints, (-delta <= u(:,i) - u(:,i-1) <= delta):'slew rate'];
end

% Construct optimizer
options = sdpsettings;
innerController = optimizer(constraints, objective, options, [x(:,1); ref; u_prev; dest], u(:,1));

% Simulate either with constant or varying reference
x0 = zeros(nx,1);
r_const = [0.8; 0.12; -0.12; pi/2];
simQuad(sys, innerController, 0, x0, T, r_const, filter, [], 1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%  Soft Constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART VII - Soft Constraints...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear previous variables.
clearvars x u objective constraints innerController ref ss xr ur dest u_prev

% Define variables for controller
ref = sdpvar(4,1);      % reference
x = sdpvar(nx,N+1);     % x
u = sdpvar(nu,N);       % u
dest = sdpvar(nx,1);    % disturbance
u_prev = sdpvar(nu,1);  % previous input

% Define system equations as x(k+1) = A x(k) + B u(k) + B_d d(k), 
% y(k) = C x(k) + C_d d(k) and the disturbance dynamics as d(k+1) = d(k).
C = eye(nx);
B_d = eye(nx);
C_d = zeros(nx);        % Set to Zero!!

% Define state observer dynamics as
% [x_hat(k+1); d_hat(k+1)] = (A_aug - L C_aug) [x_hat(k); d_hat(k)] +
% [B_aug L] [u(k) x(k)], where L is chosen such thath the state observer
% dynamics are stable and go to zero.
A_aug = [sys.A B_d; zeros(nx) eye(nx)];
B_aug = [sys.B; zeros(nx,nu)];
C_aug = [eye(nx) eye(nx)];

% Choose cost matrices for observer design
Q_ = diag([1*ones(1,nx) [50 1 1 500 10 10 0.01 ]]);
R_ = 10*eye(nx);

% For this Q_ and R_ also the offset-free trackign with r_const is feasible
% Q_ = diag(ones(2*nx,1));
% R_ = diag(ones(nx,1));
L = dlqr(A_aug',C_aug',Q_,R_)';

% Defining the filter
filter.Af = A_aug - L*C_aug;
% abs(eig(filter.Af))
filter.Bf = [B_aug L];

% Use formula of L07 slide 29 to compute xr and ur directly, which is not a
% problem, since the matrix has full rank.
Z = [sys.A-eye(nx), sys.B; C(1:4,:), zeros(nu)];
ss = Z\[-B_d*dest; ref - 0*C_d(1:4,:)*dest];
xr = ss(1:nx);
ur = ss(nx+1:end);

% Init constraints and objective function
constraints = [];
objective = 0;
for i = 1:N
    % Add state evolution constraints.
    constraints = [constraints, x(:,i+1) == sys.A*(x(:,i)) + sys.B*(u(:,i)) + B_d*dest];
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

% Define slack variable objective function parameters
S = eye(nu);        % S > 0
v = 1;              % v > lambda

% Add soft slew rate constraints and objective
constraints = [constraints, (-delta - epsilon(:,1) <= u_prev - u(:,1) <= delta + epsilon(:,1)):'soft slew rate'];
constraints = [constraints, (epsilon(:,1) >= 0):'positive slackness'];
objective = objective + v*norm(epsilon(:,1),1) + epsilon(:,1)' * S * epsilon(:,1);
for i = 2:N
% Add soft slew rate constraints
    constraints = [constraints, (-delta - epsilon(:,i) <= u(:,i) - u(:,i-1) <= delta + epsilon(:,i)):'soft slew rate'];
    constraints = [constraints, (epsilon(:,i) >= 0):'positive slackness'];
    % Add slack variables cost to objective
    objective = objective + v*norm(epsilon(:,i),1) + epsilon(:,i)' * S * epsilon(:,i);
end

% Construct optimizer
options = sdpsettings;
innerController = optimizer(constraints, objective, options, [x(:,1); ref; u_prev; dest], u(:,1));

% Simulate either with constant or varying reference
x0 = zeros(nx,1);
r_const = [0.8; 0.12; -0.12; pi/2];
[~, ~, ~, ~, deltat] = simQuad(sys, innerController, 0, x0, T, r_const, filter, [], 2);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%  FORCES Pro %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART VIII - FORCES Pro...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulate task 5 with FORCES Pro controller
x0 = zeros(nx,1);
r_const = [1.0; 0.1745; -0.1745; 1.7453];
[~, ~, ~, rt_task5_FORCES ,~] = simQuad(sys, innerController_FORCES_task5, 0, x0, T, r_const)

% Simulate task 9 with FORCES Pro controller
x0 = zeros(nx,1);
r_const = [0.8; 0.12; -0.12; pi/2];
[~, ~, ~, rt_task9_FORCES ,~] = simQuad(sys, innerController, 0, x0, T, r_const, filter, [])