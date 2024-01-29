% %% implement a car vehicle model with rear traction path tracking MPC
% % 
clearvars;
close all
clc


addpath('~/casadi-3.6.4-linux64-matlab2018b')
%addpath('/Users/mariobozza/Documents/CasADi/casadi-3.6.4-osx64-matlab2018b')
import casadi.*

%%
T = 0.1; %[s]
N = 40; % prediction horizon
lf = 0.4;
lr = 0.4;
lambda =  -1e-5;
v_max = 0.6; v_min = -v_max;
delta_max = pi/4; delta_min = -delta_max;

% States
x = SX.sym('x'); y = SX.sym('y'); psi = SX.sym('psi'); 
s = SX.sym('s');
states = [x; y; psi; s]; n_st = length(states);

% inputs
vf = SX.sym('vf'); vr = SX.sym('vr'); 
deltaf = SX.sym('deltaf'); deltar = SX.sym('deltar');
virtual_v = SX.sym('virtual_v');
controls = [vf; vr; deltaf; deltar; virtual_v]; n_con = length(controls);

%% path gen
load("x-y_path.mat")
arc_length = arc_length - arc_length(end);
%%
% kinematics
beta = atan((lf * tan(deltar) + lr * tan(deltaf)) / (lf + lr));
v = (vf * cos(deltaf) + vr * cos(deltar)) / (2 * cos(beta));

rhs = [v * cos(psi + beta); v * sin(psi + beta); ... 
       v * cos(beta) * (tan(deltaf) - tan(deltar)) / (lf + lr);
       lambda*(-52) + virtual_v]; % system r.h.s

f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)
U = SX.sym('U',n_con,N); % Decision variables (controls)
P = SX.sym('P',n_st + n_st + n_con);% Parameters (x,x_r,u_r)
X = SX.sym('X',n_st,(N+1)); % A vector that represents the states over the optimization problem.
obj = 0; % Objective function
g = [];  % constraints vector

Q = zeros(4,4); Q(1,1) = 1e7; Q(2,2) = 1e7; 
                Q(3,3) = 8e5; q(4,4) = 0.5;% weighing matrices (states)

R = zeros(5,5); R(1,1) = 1e1; R(2,2) = 1e1; ...
                R(3,3) = 1e1; R(4,4) = 1e1; ...
                R(5,5) = 0.01;% weighing matrices (controls)

W = zeros(5,5); W(1,1) = 1e8; W(2,2) = 1e8;...
                W(3,3) = 1e10; W(4,4) = 1e10;...
                W(5,5) = 1; %weight for rate of change

eps = 1e3;
st  = X(:,1); % initial state
g = [g;st-P(1:4)]; % initial condition constraints

for k = 1:N
    st = X(:,k);  con = U(:,k); 
    if k == N
        con_l = con;
    else
        con_l = U(:,k+1);   
    end
    obj = obj + (st-P(5:8))'*Q*(st-P(5:8)) + ...
          (con-P(9:13))'*R*(con-P(9:13)) + ...
          (con - con_l)'*W*(con - con_l); % calculate obj
    st_next = X(:,k+1);
    f_value = f(st,con);
    st_next_euler = st+ (T*f_value);
    g = [g;st_next-st_next_euler]; % compute constraints
end
obj = obj + eps/2*st(4)^2;

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

args.lbx(1:n_st:n_st*(N+1)-n_st+1,1) =-30; %state x lower bound
args.ubx(1:n_st:n_st*(N+1)-n_st+1,1) = 0; %state x upper bound
args.lbx(2:n_st:n_st*(N+1)-n_st+2,1) = - inf; %state y lower bound
args.ubx(2:n_st:n_st*(N+1)-n_st+2,1) = inf; %state y upper bound
args.lbx(3:n_st:n_st*(N+1)-n_st+3,1) = -inf; %state psi lower bound
args.ubx(3:n_st:n_st*(N+1)-n_st+3,1) = inf; %state psi upper bound
args.lbx(4:n_st:n_st*(N+1)-n_st+4,1) = -58; %state s lower bound
args.ubx(4:n_st:n_st*(N+1)-n_st+4,1) = 0; %state s upper bound

% input constraints
args.lbx(n_st*(N+1)+1:n_con:n_st*(N+1)+n_con*N,1) = 0;      %lb on vf
args.ubx(n_st*(N+1)+1:n_con:n_st*(N+1)+n_con*N,1) = 1;      %ub on vf
args.lbx(n_st*(N+1)+2:n_con:n_st*(N+1)+n_con*N,1) = 0;      %lb on vr
args.ubx(n_st*(N+1)+2:n_con:n_st*(N+1)+n_con*N,1) = 1;      %ub on vr
args.lbx(n_st*(N+1)+3:n_con:n_st*(N+1)+n_con*N,1) = - 0.63; %lb on deltaf
args.ubx(n_st*(N+1)+3:n_con:n_st*(N+1)+n_con*N,1) = 0.63;   %ub on deltaf
args.lbx(n_st*(N+1)+4:n_con:n_st*(N+1)+n_con*N,1) = - 0.63; %lb on deltar
args.ubx(n_st*(N+1)+4:n_con:n_st*(N+1)+n_con*N,1) = 0.63;   %ub on deltar
args.lbx(n_st*(N+1)+5:n_con:n_st*(N+1)+n_con*N,1) = 0;      %lb on virtual v
args.ubx(n_st*(N+1)+5:n_con:n_st*(N+1)+n_con*N,1) = 1;      %ub on virtual v


%----------------------------------------------
% ALL OF THE ABOVE IS JUST A PROBLEM SET UP


% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
t0 = 0;


%Initial conditions
x0 = [x_p(1); y_p(1); 0; arc_length(1)]; %initial condition state
xp0 = [x_p(1) ; y_p(1) ; 1; 0.0];    % initial condition path
xf = [x_p(end); y_p(end); 0; 0]; %last path position
s_0 = arc_length(1);
up0 = zeros(n_con,1); %initial control reference
xx(:,1) = x0; % xx contains the history of states
t(1) = t0;

u0 = zeros(N,n_con);        % two control inputs for each robot
X0 = repmat(x0,1,N+1)'; % initialization of the states decision variables

sim_tim = 100; % Maximum simulation time

% Start MPC
mpciter = 0;
xx1 = [];
u_cl=[];

% the main simulaton loop... it works as long as the error is greater
% than 10^-6 and the number of mpc steps is less than its maximum
% value.

main_loop = tic;
fin = 0; 
while(fin < 20 && mpciter < sim_tim / T )
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
    s_prev = s_0;
    s_0 = x0(4);
    if norm((x0(1:2)-xf(1:2)),2) < 5e-1
        fin = fin + 1;
    end
    % Apply the control and shift the solution
    [t0, x0, xp0, u0, up0] = PP_4ws_shift(mpciter, T, t0, x0, u, f, arc_length, s_0, s_prev, x_p, y_p);
    xx(:,mpciter+2) = x0;
    X0 = reshape(full(sol.x(1:n_st*(N+1)))',n_st,N+1)'; % get solution TRAJECTORY
    % Shift trajectory to initialize the next step
    X0 = [X0(2:end,:);X0(end,:)];
    mpciter;
    mpciter = mpciter + 1;
    error(:,mpciter) = norm((x0(1:3)-xp0(1:3)),2);
    step_time(mpciter) = toc(step);
    dim_error(:,mpciter) = sqrt((x0(1)-xp0(1))^2+(x0(2)-xp0(2))^2);
end
main_loop_time = toc(main_loop);

average_mpc_time = main_loop_time/(mpciter+1);

f = Function('f',{deltaf,deltar},{beta});
b = full(f(u_cl(:,3),u_cl(:,4)));
Draw_MPC_PP_4ws_path_trackin_carexample (t,xx,xx1,u_cl,xf,N, step_time, average_mpc_time, cost_f, dim_error, x_p, y_p, b)


