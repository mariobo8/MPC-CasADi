function [t0, x0, xp0, u0, up0] = PP_shift(mpciter, T, t0, x0, u, f, arc_length, theta_prev, theta_prev_prev, x_p, y_p)
st = x0;
con = u(1,:)';
f_value = f(st,con);
st = st+ (T*f_value);
x0 = full(st);
theta = x0(4);

x_int = interp1(arc_length, x_p, theta, 'pchip');
y_int = interp1(arc_length, y_p, theta, 'pchip');

x_int_prev = interp1(arc_length, x_p, theta_prev, 'pchip');
y_int_prev = interp1(arc_length, y_p, theta_prev, 'pchip');

x_int_prev_prev = interp1(arc_length, x_p, theta_prev_prev, 'pchip');
y_int_prev_prev = interp1(arc_length, y_p, theta_prev_prev, 'pchip');
if mpciter == 0
    psi_int = 0;
    psi_int_prev = 0;
else
    psi_int = atan2((y_int - y_int_prev) , (x_int - x_int_prev));
    psi_int_prev = atan2((y_int_prev - y_int_prev_prev), ...
         (x_int_prev - x_int_prev_prev));
end

%compure ref vector

xp0 = [x_int; y_int; psi_int; 0];
t0 = t0 + T;
u0 = [u(2:size(u,1),:);u(size(u,1),:)];


%this above is a super big approximation

up0 = [((x_int-x_int_prev)-(y_int-y_int_prev))/(T*(cos(x0(3))+sin(x0(3)))); 
        atan2((psi_int - psi_int_prev),(T*((x_int-x_int_prev)-(y_int-y_int_prev))...
        /(T*(cos(x0(3))+sin(x0(3))))));
        0];
% %up0 will be the vector of uref
% subplot(3,1,1)
% scatter(t0,x0(1))
% hold on
% scatter(t0,xp0(1))
% hold on
% subplot(3,1,2)
% scatter(t0,x0(2))
% hold on
% scatter(t0,xp0(2))
% hold on
% subplot(3,1,3)
% scatter(t0,x0(3))
% hold on
% scatter(t0,xp0(3))
% hold on
end