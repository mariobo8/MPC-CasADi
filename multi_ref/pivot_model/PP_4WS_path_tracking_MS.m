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
N = 10; % prediction horizon
lf = 0.4;
lr = 0.4;
lambda =  -1e-1;
v_max = 0.6; v_min = -v_max;
delta_max = pi/4; delta_min = -delta_max;

% States
x = SX.sym('x'); y = SX.sym('y'); psi = SX.sym('psi'); 
s = SX.sym('s');
states = [x; y; psi; s]; n_st = length(states);

% inputs
vf = SX.sym('vf'); vr = SX.sym('vr'); 
deltaf = SX.sym('deltaf'); deltar = SX.sym('deltar');
alpha = SX.sym('alpha');
virtual_v = SX.sym('virtual_v');
controls = [vf; vr; deltaf; deltar; alpha ;virtual_v]; n_con = length(controls);

%% path gen
load("s_shape_path.mat")
% x_p = x_p(1:10);
% y_p = x_p(1:10);
% arc_length = arc_length(1:10);
arc_length = arc_length;% - arc_length(end);
%%
% kinematics
beta = atan2((lf * tan(deltar) + lr * (tan(deltaf)*cos(alpha) + sin(alpha))) , (lf + lr*(cos(alpha) - tan(deltaf)*sin(alpha))));
v = (vf * cos(deltaf + alpha) + vr * cos(deltar)) / (2 * cos(beta));
a = 3;

rhs = [v * cos(psi + beta); v * sin(psi + beta); ... 
       v * (tan(deltaf) * cos(alpha + beta) + sin(beta) * ...
       (1 + cos(alpha)) - tan(deltar) * cos(beta)) / (lf + lr);
       virtual_v]; % system r.h.s

f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)
U = SX.sym('U',n_con,N); % Decision variables (controls)
P = SX.sym('P',n_st + n_st*N + n_con*N);% Parameters (x,x_r,u_r)
X = SX.sym('X',n_st,(N+1)); % A vector that represents the states over the optimization problem.
obj = 0; % Objective function
g = [];  % constraints vector

Q = zeros(4,4); Q(1,1) = 4e5; Q(2,2) = 4e5; 
                Q(3,3) = 4e4; Q(4,4) = 0;% weighing matrices (states)

R = zeros(6,6); R(1,1) = 1e2; R(2,2) = 1e2; ...
                R(3,3) = 1e2; R(4,4) = 1e2; ...
                R(5,5) = 1e2; R(6,6) = 0;% weighing matrices (controls)

% W = zeros(5,5); W(1,1) = 1e7; W(2,2) = 1e7;...
%                 W(3,3) = 1e4; W(4,4) = 1e4;...
%                 W(5,5) = 1e5; %weight for rate of change


W = zeros(6,6); W(1,1) = 1e5; W(2,2) = 1e5;...
                W(3,3) = 4e4; W(4,4) = 4e4;...
                W(5,5) = 4e4; W(6,6) = 1e1; %weight for rate of change
eps = 1e3;
st  = X(:,1); % initial state
g = [g;st-P(1:4)]; % initial condition constraints

for k = 1:N
    st = X(:,k);  con = U(:,k); conn = [U(5,k); U(5,k); 0; 0; 0; 0];
    ind_st = (((k-1)*4+5):((k-1)*4+8));
    ind_con = ((N-1)*4+9+(k-1)*n_con : (N-1)*4+9+(k-1)*n_con + 5);
    obj = obj + (st-P(ind_st))'*Q*(st-P(ind_st)) + ...
          (con)'*R*(con) + ...
          (con - P(ind_con))'*W*(con - P(ind_con)); % calculate obj
          %([U(5,k); U(5,k)] - [U(1,k); U(2,k)])'*F*([U(5,k); U(5,k)] - [U(1,k); U(2,k)]) +
          
    st_next = X(:,k+1);
    f_value = f(st,con);
    st_next_euler = st+ (T*f_value);
    g = [g;st_next-st_next_euler]; % compute constraints
end
 obj = obj - eps/2*(st(4))^2;

% make the decision variable one column  vector
OPT_variables = [reshape(X,n_st*(N+1),1);reshape(U,n_con*N,1)];

nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 6000;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;



solver = nlpsol('solver', 'ipopt', nlp_prob, opts);

args = struct;

args.lbg(1:n_st*(N+1)) = 0;  % -1e-20  % Equality constraints
args.ubg(1:n_st*(N+1)) = 0;  % 1e-20   % Equality constraints

args.lbx(1:n_st:n_st*(N+1),1) =-inf; %state x lower bound
args.ubx(1:n_st:n_st*(N+1),1) = 0; %state x upper bound
args.lbx(2:n_st:n_st*(N+1),1) = - inf; %state y lower bound
args.ubx(2:n_st:n_st*(N+1),1) = inf; %state y upper bound
args.lbx(3:n_st:n_st*(N+1),1) = -inf; %state psi lower bound
args.ubx(3:n_st:n_st*(N+1),1) = inf; %state psi upper bound
args.lbx(4:n_st:n_st*(N+1),1) = -inf; %state s lower bound
args.ubx(4:n_st:n_st*(N+1),1) = arc_length(end); %state s upper bound

% input constraints
args.lbx(n_st*(N+1)+1:n_con:n_st*(N+1)+n_con*N,1) = -1;      %lb on vf
args.ubx(n_st*(N+1)+1:n_con:n_st*(N+1)+n_con*N,1) = 1;      %ub on vf
args.lbx(n_st*(N+1)+2:n_con:n_st*(N+1)+n_con*N,1) = -1;      %lb on vr
args.ubx(n_st*(N+1)+2:n_con:n_st*(N+1)+n_con*N,1) = 1;      %ub on vr
args.lbx(n_st*(N+1)+3:n_con:n_st*(N+1)+n_con*N,1) = - 1.05; %lb on deltaf
args.ubx(n_st*(N+1)+3:n_con:n_st*(N+1)+n_con*N,1) = 1.05;   %ub on deltaf
args.lbx(n_st*(N+1)+4:n_con:n_st*(N+1)+n_con*N,1) = - 1.05 ; %lb on deltar
args.ubx(n_st*(N+1)+4:n_con:n_st*(N+1)+n_con*N,1) = 1.05;   %ub on deltar
args.lbx(n_st*(N+1)+5:n_con:n_st*(N+1)+n_con*N,1) = -1.05;      %lb on virtual alpha
args.ubx(n_st*(N+1)+5:n_con:n_st*(N+1)+n_con*N,1) = 1.05;      %ub on virtual alpha
args.lbx(n_st*(N+1)+6:n_con:n_st*(N+1)+n_con*N,1) = 0;      %lb on virtual v
args.ubx(n_st*(N+1)+6:n_con:n_st*(N+1)+n_con*N,1) = 1;      %ub on virtual v


%----------------------------------------------
% ALL OF THE ABOVE IS JUST A PROBLEM SET UP


% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
t0 = 0;


%Initial conditions
x0 = [x_p(1); y_p(1); 0; arc_length(1)]; %initial condition state
xp0 = [];
for jj = 1 : N  
    xp0 = [xp0; x_p(1) ; y_p(1) ; 0.0; 0.0];    % initial condition path
end
xf = [x_p(end); y_p(end); 0; 0]; %last path position
s_0 = arc_length(1);
up0 = zeros(n_con*N,1); %initial control reference
xx(:,1) = x0; % xx contains the history of states
t(1) = t0;

u0 = zeros(N,n_con); 
X0 = repmat(x0,1,N+1)'; % initialization of the states decision variables

sim_tim = 15; % Maximum simulation time

% Start MPC
mpciter = 0;
xx1 = [];
u_cl=[];


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
    x_s = xx1(end,1,mpciter+1);
    u_cl= [u_cl ; u(1,:)];
    t(mpciter+1) = t0;
    s_prev = s_0;
    s_0 = x0(4);
    if norm((x0(1:2)-xf(1:2)),2) < 5e-1
        fin = fin + 1;
    end
    s_00 = xx1(end,4,end);
    % Apply the control and shift the solution
    [t0, x0, xp0, u0, up0, ref] = PP_4ws_shift(mpciter, T, t0, x0, u, f, arc_length, s_0, x_p, y_p, N, x_s, x0(1),x0(2),s_00, n_con);
    x_r(:,:,mpciter+1) = ref;
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

f = Function('f',{deltaf,deltar,alpha},{beta});
b = full(f(u_cl(:,3),u_cl(:,4),u_cl(:,5)));
Draw_MPC_PP_4ws_path_trackin_carexample (t,xx,xx1,u_cl,xf,N, step_time, average_mpc_time, cost_f, dim_error, x_p, y_p, b, arc_length, x_r, a)

%%
omega_f = diff(u_cl(:,3));
omega_r = diff(u_cl(:,4));
omega_p = diff(u_cl(:,5));
energy_f = sum(abs(omega_f)*T);
energy_r = sum(abs(omega_r)*T);
energy_p = sum(abs(omega_p)*T);
tot_energy=energy_r+energy_f+energy_p;
figure
plot(t(1:end-1),omega_f)
hold on
plot(t(1:end-1), omega_r)
hold on 
plot(t(1:end-1), omega_p)
grid on
disp(energy_f)
disp(energy_r)
disp(energy_p)
disp(tot_energy)
