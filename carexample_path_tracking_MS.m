%% implement a car vehicle model with rear traction path tracking MPC
% 
% clearvars;
% close all
% clc


addpath('~/casadi-3.6.4-linux64-matlab2018b')

import casadi.*

%%
T = 0.3; %[s]
N = 30; % prediction horizon
lf = 0.4;
lr = 0.4;
lambda = - 1e-3;
v_max = 0.6; v_min = -v_max;
delta_max = pi/4; delta_min = -delta_max;

% States
x = SX.sym('x'); y = SX.sym('y'); psi = SX.sym('psi'); 
theta = SX.sym('theta');
states = [x; y; psi; theta]; n_st = length(states);

% inputs
v = SX.sym('v'); 
delta = SX.sym('delta'); 
virtual_v = SX.sym('virtual_v');
controls = [v; delta; virtual_v]; n_con = length(controls);

%% path gen
mu = [6; 20; 5; 0.35];  
rho = - mu(1) * log(mu(2)./(mu(3) + abs(theta))) .* sin(mu(4) * theta);

%%
% kinematics

rhs = [v * cos(psi); v * sin(psi); ... 
       v * tan(delta);
       -lambda * theta + virtual_v]; % system r.h.s

f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)
U = SX.sym('U',n_con,N); % Decision variables (controls)
P = SX.sym('P',n_st + n_st + n_con);% Parameters (x,x_r,u_r)
X = SX.sym('X',n_st,(N+1)); % A vector that represents the states over the optimization problem.
obj = 0; % Objective function
g = [];  % constraints vector

% Q = zeros(4,4); Q(1,1) = 8e6; Q(2,2) = 2e6; 
%                 Q(3,3) = 8e5; q(4,4) = 10;% weighing matrices (states)
% R = zeros(3,3); R(1,1) = 10; R(2,2) = 1; ...
%                 R(3,3) = 1;% weighing matrices (controls)
% eps = 1500;
st  = X(:,1); % initial state
g = [g;st-P(1:4)]; % initial condition constraints

for mu = 1:N
    st = X(:,mu);  con = U(:,mu); 
    obj = obj + (st-P(5:8))'*Q*(st-P(5:8)) + (con-P(9:11))'*R*(con-P(9:11)) + eps/2*st(4)^2; % calculate obj
    st_next = X(:,mu+1);
    f_value = f(st,con);
    st_next_euler = st+ (T*f_value);
    g = [g;st_next-st_next_euler]; % compute constraints
end
% make the decision variable one column  vector
OPT_variables = [reshape(X,n_st*(N+1),1);reshape(U,n_con*N,1)];

nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 2000;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;



solver = nlpsol('solver', 'ipopt', nlp_prob, opts);

args = struct;

args.lbg(1:n_st*(N+1)) = 0;  % -1e-20  % Equality constraints
args.ubg(1:n_st*(N+1)) = 0;  % 1e-20   % Equality constraints

%parameters for constr
eps = 0; %also employed to construct an extendedset 
            % Z_ε (introduced later) taking the particular structure of
            % the path-following problem into account.

%
args.lbx(1:n_st:n_st*(N+1)-n_st+1,1) =-30; %state x lower bound
args.ubx(1:n_st:n_st*(N+1)-n_st+1,1) = -eps; %state x upper bound
args.lbx(2:n_st:n_st*(N+1)-n_st+2,1) = - inf; %state y lower bound
args.ubx(2:n_st:n_st*(N+1)-n_st+2,1) = inf; %state y upper bound
args.lbx(3:n_st:n_st*(N+1)-n_st+3,1) = -inf; %state psi lower bound
args.ubx(3:n_st:n_st*(N+1)-n_st+3,1) = inf; %state psi upper bound
args.lbx(4:n_st:n_st*(N+1)-n_st+4,1) = -30; %state theta lower bound
args.ubx(4:n_st:n_st*(N+1)-n_st+4,1) = 0; %state theta upper bound

% input constraints
args.lbx(n_st*(N+1)+1:n_con:n_st*(N+1)+n_con*N-n_con+1,1) = 0; %lb on v
args.ubx(n_st*(N+1)+1:n_con:n_st*(N+1)+n_con*N-n_con+1,1) = 1; %ub on v
args.lbx(n_st*(N+1)+2:n_con:n_st*(N+1)+n_con*N-n_con+2,1) = - 0.63; %lb on delta
args.ubx(n_st*(N+1)+2:n_con:n_st*(N+1)+n_con*N-n_con+2,1) = 0.63; %ub on delta
args.lbx(n_st*(N+1)+3:n_con:n_st*(N+1)+n_con*N-n_con+3,1) = 0; %lb on virtual v
args.ubx(n_st*(N+1)+3:n_con:n_st*(N+1)+n_con*N-n_con+3,1) = 1; %ub on virtual v


%----------------------------------------------
% ALL OF THE ABOVE IS JUST A PROBLEM SET UP


% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
t0 = 0;

%function to compute values of the path
f2 = Function('f',{theta},{rho});
rho_0 = full(f2(-30));
rho_f = full(f2(0));
%Initial conditions
x0 = [-30; 3; 0.0; -30]; %initial condition state
xp0 = [-30 ; rho_0 ; -1; 0.0];    % initial condition path
xf = [0; rho_f; 0; 0]; %last path position

up0 = [0; 0; 0]; %initial control reference
xx(:,1) = x0; % xx contains the history of states
t(1) = t0;

u0 = zeros(N,n_con);        % two control inputs for each robot
X0 = repmat(x0,1,N+1)'; % initialization of the states decision variables

sim_tim = 15; % Maximum simulation time

% Start MPC
mpciter = 0;
xx1 = [];
u_cl=[];

% the main simulaton loop... it works as long as the error is greater
% than 10^-6 and the number of mpc steps is less than its maximum
% value.
theta_prev = -30;
main_loop = tic;
while(norm((x0-xf),2) > 1e-1 && mpciter < sim_tim / T)
    step = tic;
    args.p   = [x0;xp0;up0]; % set the values of the parameters vector
    % initial value of the optimization variables
    args.x0  = [reshape(X0',n_st*(N+1),1);reshape(u0',n_con*N,1)];
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
    cost_f(mpciter+1) = full(sol.f);
    u = reshape(full(sol.x(n_st*(N+1)+1:end))',n_con,N)'; % get controls only from the solution
    xx1(:,1:n_st,mpciter+1)= reshape(full(sol.x(1:n_st*(N+1)))',n_st,N+1)'; % get solution TRAJECTORY
    u_cl= [u_cl ; u(1,:)];
    t(mpciter+1) = t0;
    theta_prev_prev = theta_prev;
    theta_prev = x0(4);
   
    % Apply the control and shift the solution
    [t0, x0, xp0, u0, up0] = shift(T, t0, x0, u, f, f2, theta_prev, theta_prev_prev);
    xx(:,mpciter+2) = x0;
    X0 = reshape(full(sol.x(1:n_st*(N+1)))',n_st,N+1)'; % get solution TRAJECTORY
    % Shift trajectory to initialize the next step
    X0 = [X0(2:end,:);X0(end,:)];
    mpciter;
    mpciter = mpciter + 1;
    error(:,mpciter) = norm((x0(1:3)-xp0(1:3)),2);
    step_time(mpciter) = toc(step);
end
main_loop_time = toc(main_loop);

average_mpc_time = main_loop_time/(mpciter+1);


%Draw_MPC_path_trackin_carexample (t,xx,xx1,u_cl,xf,N, step_time, average_mpc_time, f2, cost_f)


