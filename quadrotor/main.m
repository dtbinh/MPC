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

%%%%%%%%%%%%%%%%%%%%%    First MPC controller %%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART I - First MPC controller...\n')

%Define an instance of the class LTIsystem using the give parameters
system = LTISystem('A',sys.A,'B',sys.B,'Ts',sys.Ts);

Q = eye(size(sys.A));
R = eye(size(sys.B,2));
N = 10;
%Get the LQR-Control law and the infinit-horizon-solution to the discrete
%ARE
[K,P_infty,~] = dlqr(sys.A,sys.B,Q,R);

%Define the LQR-Closed-Loop system and use the memeberfunction
%invariantSet() to find the maximum LQR-Control-Invariant Set which will be
%used as terminal constraint.
systemC = LTISystem('A',sys.A-sys.B*K,'Ts',sys.Ts);
Xf = systemC.invariantSet('maxIterations',30);
%Take the input-constraints of the inner loop into account!
A_inputConstrains = [...
    1 0 0 0 0 0 0;...
    -1 0 0 0 0 0 0;...
     0 1 0 0 0 0 0;...
     0 -1 0 0 0 0 0;...
     0 0 1 0 0 0 0;...
     0 0 -1 0 0 0 0;...
     0 0 0 0 1 0 0;...
     0 0 0 0 -1 0 0;...
     0 0 0 0 0 1 0;...
     0 0 0 0 0 -1 0;...
     0 0 0 0 0 0 1;...
     0 0 0 0 0 0 -1];
 b_inputConstraints = [...
     1;...
     1;...
     10/180*pi;...
     10/180*pi;...
     10/180*pi;...
     10/180*pi;...
     15/180*pi;...
     15/180*pi;...
     15/180*pi;...
     15/180*pi;...
     60/180*pi;...
     60/180*pi];
IC = polytope(A_inputConstrains,b_inputConstraints);   

%Defining the necessary sdp-variables for all N iterations
nx = size(sys.A,1);
nu = size(sys.B,2);
x = sdpvar(nx,1);
u = sdpvar(repmat(nu,1,N),ones(1,N));
%Defining the constraints and the objective for all iterations
constraints = [];
objective = 0;
for k = 1:N-1
    %Equality constraint
    x = sys.A*x + sys.B*u{k};
    %Update the objective
    objective = objective + x'*Q*x + u{k}'*R*u{k};
    %Update the constraints: Inner loop state constraints
    constraints = [constraints,...
        abs(x(1))<=1,...
        abs(x(2))<=10/180*pi,...
        abs(x(3))<=10/180*pi,...
        abs(x(5))<=15/180*pi,...
        abs(x(6))<=15/180*pi,...
        abs(x(7))<=60/180*pi];
    %Update the constraints: Inner loop input constraints
    constraints = [constraints,...
        u{k}>zeros(4,1),...
        u{k}<ones(4,1)]
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
