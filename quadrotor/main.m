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
x0 = [-1; 0.1745; -0.1745; 0.8727; 0; 0; 0];    % Initial condition

%%%%%%%%%%%%%%%%%%%%%    First MPC controller %%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART I - First MPC controller...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Chose parameters N, Q, R, P. Chose final cost to be the the solution of
% DARE of the unconstrained LQR control law.
N = 20;
Q = 10*eye(nx);
R = eye(nu);
[K_inf, P, ~] = dlqr(sys.A, sys.B, Q, R);

% Compute the maximum invariant set X_f for the cl-system 
% x_k+1 = (A - B*K_inf)*x_k.
clSystem = LTISystem('A', sys.A - sys.B*K_inf, 'Ts', sys.Ts);
X_cl = clSystem.invariantSet();
A_X_f = X_cl.A;
b_X_f = X_cl.b;

% Define constraints on states and inputs.
stateConstraint = [1; 10/180*pi; 10/180*pi; Inf; 15/180*pi; 15/180*pi; 60/180*pi];

% Define optimization variables x and u in appropriate sizes.
x = sdpvar(nx,N+1);   % [dzdt; alpha; beta; gamma; dalphadt, dbetadt, dgammadt]
u = sdpvar(nu,N);   % [u_1; u_2; u_3; u_4]

% Generate constraints and objective function.
constraints = [];
objective = 0;
for i = 1:N
    % Add state evolution constraints.
    constraints = [constraints, x(:,i+1) == sys.A*x(:,i) + sys.B*u(:,i)];
    % Add state constraints.
    constraints = [constraints, -stateConstraint <= x(:,i) <= stateConstraint];
    % Add input constraints.
    constraints = [constraints, zeros(nu,1) <= u(:,i) + us <= ones(nu,1)];
    
    % Add to objective function
    objective = objective + x(:,i)' * Q * x(:,i) + u(:,i)' * R * u(:,i);
end
% Add final state constraints and objective function.
constraints = [constraints, A_X_f * x(:,N+1) <= b_X_f];
objective = x(:,N+1)' * P * x(:,N+1);

% Call the optimizer
options = sdpsettings;
innerController = optimizer(constraints, objective, options, x(:,1), u(:,1));
simQuad(sys, innerController, bForces, x0, T);

%%%%%%%%%%%%%%%%%%%%%  Reference Tracking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART II - Reference tracking...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%  First simulation of the nonlinear model %%%%%%%%%%%%%%%%%
fprintf('PART III - First simulation of the nonlinear model...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%  Offset free MPC  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
