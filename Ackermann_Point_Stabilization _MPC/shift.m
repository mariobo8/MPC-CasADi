function [t0, x0, xp0, u0, up0] = shift(T, t0, x0, u, f, f2, theta_prev, theta_prev_prev)
st = x0;
con = u(1,:)';
f_value = f(st,con);
st = st+ (T*f_value);
x0 = full(st);
theta = x0(4);
rho = full(f2(theta));
rho_prev = full(f2(theta_prev));
rho_prev_prev = full(f2(theta_prev));
%compure ref vector
drho_dtheta = (rho - rho_prev)/...
              (theta - theta_prev);

d2rho_dtheta2 = (rho - 2 * (rho - rho_prev) + (rho - rho_prev_prev))/...
                ((theta - theta_prev)*(theta_prev-theta_prev_prev));
xp0 = [theta; rho; atan(drho_dtheta); 0];
t0 = t0 + T;
%add last one twice to make them the same size?!
u0 = [u(2:size(u,1),:);u(size(u,1),:)];


%this above is a super big approximation

up0 = [(theta-theta_prev)/T*sqrt(1+(drho_dtheta)^2); 
        atan((1+(drho_dtheta^2)^(-3/2))*d2rho_dtheta2);
        0];
%up0 will be the vector of uref
end