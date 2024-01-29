% point stabilization Ackermann + Single shooting
clear all
close all
clc
addpath('/Users/mariobozza/Documents/CasADi/casadi-3.6.4-osx64-matlab2018b')

import casadi.*

T = 0.2; % sampling time [s]
N = 40; % prediction horizon
rob_diam = 0.3;
lf = 0.4;
lr = 0.4;
v_max = 0.6; v_min = -v_max;
delta_max = pi/4; delta_min = -delta_max;


% States
x = SX.sym('x'); y = SX.sym('y'); psi = SX.sym('psi');
states = [x; y; psi]; n_states = length(states);

% inputs
vf = SX.sym('vf'); vr = SX.sym('vr'); 
deltaf = SX.sym('deltaf'); deltar = SX.sym('deltar');

controls = [vf; vr; deltaf; deltar]; n_controls = length(controls);

% kinematics
beta = atan((lf * tan(deltar) + lr * tan(deltaf)) / (lf + lr));
v = (vf * cos(deltaf) + vr * cos(deltar)) / (2 * cos(beta));

rhs = [v * cos(psi + beta); v * sin(psi + beta); ... 
       v * cos(beta) * (tan(deltaf) + tan(deltar)) / (lf + lr)]; % system r.h.s

f = Function('f', {states,controls}, {rhs}); % nonlinear mapping function f(x,u)
U = SX.sym('U', n_controls, N); % Decision variables (controls)
P = SX.sym('P', n_states + n_states);
% parameters (which include the initial and the reference state of the robot)

X = SX.sym('X', n_states, (N+1));
% A Matrix that represents the states over the optimization problem.
%% thist part is skipped in multiple shooting
% compute solution symbolically
X(:,1) = P(1:n_states); % initial state
for k = 1:N
    st = X(:,k);  con = U(:,k);
    f_value  = f(st,con);
    st_next  = st+ (T*f_value);
    X(:,k+1) = st_next;
end
% this function to get the optimal trajectory knowing the optimal solution
ff=Function('ff',{U,P},{X});
%%
obj = 0; % Objective function
g = [];  % constraints vector

Q = zeros(3,3); Q(1,1) = 1;Q(2,2) = 5;Q(3,3) = 0.1; % weighing matrices (states)
R = zeros(4,4); R(1,1) = 0.5; R(2,2) = 0.5; ...
                R(3,3) = 0.1; R(4,4) = 0.1;% weighing matrices (controls)

% compute objective
for k=1:N
    st = X(:,k);  con = U(:,k);
    obj = obj+(st-P(4:6))'*Q*(st-P(4:6)) + con'*R*con; % calculate obj
end

% compute constraints
for k = 1:N+1   % box constraints due to the map margins
    g = [g ; X(1,k)];   %state x
    g = [g ; X(2,k)];   %state y
end

% make the decision variables one column vector
OPT_variables = reshape(U,n_controls*N,1);
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);


args = struct;
% inequality constraints (state constraints)
args.lbg = -2;  % lower bound of the states x and y
args.ubg = 2;   % upper bound of the states x and y 

% input constraints
args.lbx(1:4:n_controls*N-4,1) = v_min; args.lbx(2:4:n_controls*N-3,1)   = v_min;
args.ubx(1:4:n_controls*N-4,1) = v_max; args.ubx(2:4:n_controls*N-3,1)   = v_max;
args.lbx(3:4:n_controls*N-1,1) = delta_min; args.lbx(4:4:n_controls*N,1)   = delta_min;
args.ubx(3:4:n_controls*N-1,1) = delta_max; args.ubx(4:4:n_controls*N,1)   = delta_max;


%----------------------------------------------
% ALL OF THE ABOVE IS JUST A PROBLEM SETTING UP


% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
t0 = 0;
x0 = [0 ; 0 ; 0.0];    % initial condition.
xs = [1.5 ; 1.5 ; 0]; % Reference posture.

xx(:,1) = x0; % xx contains the history of states
t(1) = t0;

u0 = zeros(N,4);  % four control inputs 

sim_tim = 20; % Maximum simulation time

% Start MPC
mpciter = 0;
xx1 = [];
u_cl=[];


% the main simulaton loop... it works as long as the error is greater
% than 10^-2 and the number of mpc steps is less than its maximum
% value.
main_loop = tic;
while(norm((x0-xs),2) > 1e-2 && mpciter < sim_tim / T)
    args.p   = [x0;xs]; % set the values of the parameters vector
    args.x0 = reshape(u0',n_controls*N,1); % initial value of the optimization variables
    %tic
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
            'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);    
    %toc
    u = reshape(full(sol.x)',n_controls,N)';
    ff_value = ff(u',args.p); % compute OPTIMAL solution TRAJECTORY
    xx1(:,1:3,mpciter+1)= full(ff_value)';
    
    u_cl= [u_cl ; u(1,:)];
    t(mpciter+1) = t0;
    [t0, x0, u0] = shift(T, t0, x0, u,f); % get the initialization of the next optimization step
    
    xx(:,mpciter+2) = x0;  
    mpciter
    mpciter = mpciter + 1;
end;
main_loop_time = toc(main_loop);
ss_error = norm((x0-xs),2)
average_mpc_time = main_loop_time/(mpciter+1)

Draw_MPC_point_stabilization_4WS (t,xx,xx1,u_cl,xs,N,rob_diam) % a drawing function


