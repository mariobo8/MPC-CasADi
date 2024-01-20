% first casadi test for mpc fpr mobile robots
clear all
close all
clc

addpath('~/casadi-3.6.4-linux64-matlab2018b')
import casadi.*

T = 0.2; %[s]
N = 30; % prediction horizon
rob_diam = 0.3;
lf = 0.4;
lr = 0.4;
m = 0.4; % boogie halfh length
v_max = 0.5; v_min = -v_max;
delta_max = pi/4; delta_min = -delta_max;
alpha_max = pi/6; alpha_min = -alpha_max; 

% States
x = SX.sym('x'); y = SX.sym('y'); psi = SX.sym('psi');
states = [x; y; psi]; n_states = length(states);

% inputs
vf = SX.sym('vf'); vr = SX.sym('vr'); 
deltaf = SX.sym('deltaf'); deltar = SX.sym('deltar');
alpha = SX.sym('alpha');

controls = [vf; vr; deltaf; deltar; alpha]; n_controls = length(controls);

% kinematics
beta = atan(((lf - m * sin(alpha))* tan(deltar) + lr * tan(deltaf + alpha)) / ((lf - m * sin(alpha)) + lr));
v = (vf * cos(deltaf + alpha) + vr * cos(deltar)) / (2 * cos(beta));

rhs = [v * cos(psi + beta); v * sin(psi + beta); ... 
       v * cos(beta) * (tan(deltaf + alpha) - tan(deltar)) / ((lf - m * sin(alpha)) + lr)]; % system r.h.s

f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)
U = SX.sym('U',n_controls,N); % Decision variables (controls)
P = SX.sym('P',n_states + n_states);
% parameters (which include the initial state and the reference state)

X = SX.sym('X',n_states,(N+1));
% A vector that represents the states over the optimization problem.

obj = 0; % Objective function
g = [];  % constraints vector

Q = zeros(3,3); Q(1,1) = 15;Q(2,2) = 15;Q(3,3) = 1; % weighing matrices (states)

R = zeros(5,5); R(1,1) = 2.5; R(2,2) = 2.5; ...
                R(3,3) = 0.1; R(4,4) = 0.1;
                R(5,5) = 0.1;% weighing matrices (controls)

W = 1; % wieghting matrix for side sliding
G = zeros(5,5); G(1,1) = 4; G(2,2) = 4; 
                G(3,3) = 4; G(3,3) = 4;% weighing matrices (acceleration)
                G(5,5) = 4;
st  = X(:,1); % initial state
g = [g;st-P(1:3)]; % initial condition constraints

for k = 1:N

    st = X(:,k);  con = U(:,k); 
    if k == N
        con_l = con;
    elseif k == 1
        con_l = U(:,k);
    else
        con_l = U(:,k-1); 
    end
    obj = obj+(st-P(4:6))'*Q*(st-P(4:6)) + con'*R*con ...
        + (con(4)+con(3))'*W*(con(4)+con(3)) ... %to be fixed
        + (con-con_l)'*G*(con-con_l) ; % calculate obj
    st_next = X(:,k+1);
    f_value = f(st,con);
    st_next_euler = st+ (T*f_value);
    g = [g;st_next-st_next_euler]; % compute constraints
end
% make the decision variable one column  vector
OPT_variables = [reshape(X,3*(N+1),1);reshape(U,n_controls*N,1)];

nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 2000;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

args = struct;

args.lbg(1:3*(N+1)) = 0;  % -1e-20  % Equality constraints
args.ubg(1:3*(N+1)) = 0;  % 1e-20   % Equality constraints

args.lbx(1:3:3*(N+1),1) = -inf; %state x lower bound
args.ubx(1:3:3*(N+1),1) = inf; %state x upper bound
args.lbx(2:3:3*(N+1),1) = -inf; %state y lower bound
args.ubx(2:3:3*(N+1),1) = inf; %state y upper bound
args.lbx(3:3:3*(N+1),1) = -inf; %state theta lower bound
args.ubx(3:3:3*(N+1),1) = inf; %state theta upper bound

% input constraints
args.lbx(3*(N+1)+1:n_controls:3*(N+1)+n_controls*N-4,1) = v_min; 
args.lbx(3*(N+1)+2:n_controls:3*(N+1)+n_controls*N-3,1)   = v_min;
args.ubx(3*(N+1)+1:n_controls:3*(N+1)+n_controls*N-4,1) = v_max;
args.ubx(3*(N+1)+2:n_controls:3*(N+1)+n_controls*N-3,1)   = v_max;
args.lbx(3*(N+1)+3:n_controls:3*(N+1)+n_controls*N-2,1) = delta_min; 
args.lbx(3*(N+1)+4:n_controls:3*(N+1)+n_controls*N-1,1)   = delta_min;
args.ubx(3*(N+1)+3:n_controls:3*(N+1)+n_controls*N-2,1) = delta_max;
args.ubx(3*(N+1)+4:n_controls:3*(N+1)+n_controls*N-1,1)   = delta_max;
%% const on pivot
args.lbx(3*(N+1)+5:n_controls:3*(N+1)+n_controls*N,1) = alpha_min;
args.ubx(3*(N+1)+5:n_controls:3*(N+1)+n_controls*N,1)   = alpha_max;

%----------------------------------------------
% ALL OF THE ABOVE IS JUST A PROBLEM SET UP


% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
t0 = 0;
x0 = [0 ; 0 ; 0.0];    % initial condition.
xs = [2 ; 2 ; pi/2]; % Reference posture.

xx(:,1) = x0; % xx contains the history of states
t(1) = t0;

u0 = zeros(N,n_controls);        % two control inputs for each robot
X0 = repmat(x0,1,N+1)'; % initialization of the states decision variables

sim_tim = 10; % Maximum simulation time

% Start MPC
mpciter = 0;
xx1 = [];
u_cl=[];

% the main simulaton loop... it works as long as the error is greater
% than 10^-6 and the number of mpc steps is less than its maximum
% value.
main_loop = tic;
while(norm((x0-xs),2) > 1e-2 && mpciter < sim_tim / T)
    args.p   = [x0;xs]; % set the values of the parameters vector
    % initial value of the optimization variables
    args.x0  = [reshape(X0',3*(N+1),1);reshape(u0',n_controls*N,1)];
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);

    u = reshape(full(sol.x(3*(N+1)+1:end))',n_controls,N)'; % get controls only from the solution
    xx1(:,1:3,mpciter+1)= reshape(full(sol.x(1:3*(N+1)))',3,N+1)'; % get solution TRAJECTORY
    u_cl= [u_cl ; u(1,:)];
    t(mpciter+1) = t0;
    % Apply the control and shift the solution
    [t0, x0, u0] = shift(T, t0, x0, u,f);
    xx(:,mpciter+2) = x0;
    X0 = reshape(full(sol.x(1:3*(N+1)))',3,N+1)'; % get solution TRAJECTORY
    % Shift trajectory to initialize the next step
    X0 = [X0(2:end,:);X0(end,:)];
    mpciter;
    mpciter = mpciter + 1;
end
main_loop_time = toc(main_loop);
ss_error = norm((x0-xs),2)
average_mpc_time = main_loop_time/(mpciter+1)

f = Function('f',{deltaf,deltar,alpha},{beta});
b = full(f(u_cl(:,3),u_cl(:,4),u_cl(:,5))); %compute value for beta

Draw_MPC_point_stabilization_4WS (t,xx,xx1,u_cl,xs,N,rob_diam, b)


