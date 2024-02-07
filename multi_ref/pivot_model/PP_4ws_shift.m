function [t0, x0, xp0, u0, up0, ref] = PP_4ws_shift(mpciter, T, t0, x0, u, f, arc_length, s_prev, x_p, y_p, N, x_s, x0_1,y0_1,s_00, n_con)
st = x0;
con = u(1,:)';
f_value = f(st,con);
st = st+ (T*f_value);
x0 = full(st);
s = x0(4);

    index = dsearchn([x_p,y_p],[x0_1,y0_1]);
    s_0 = arc_length(index); %from the actual position
    % scatter(mpciter, s_0, 'k')
    % hold on
    % scatter(mpciter, s, 'g')
    % hold on
    % legend('actual pos','goal')
    delta_s = (s - s_0)/N;
    sjj = s_0;
    s_prev = s_0;
    ref = [];
    xp0=[];
for jj=1:N-1
    sjj = s_0 + (jj) * delta_s;
    x_int = interp1(arc_length, x_p, sjj, 'pchip');
    y_int = interp1(arc_length, y_p, sjj, 'pchip');
    
    x_int_prev = interp1(arc_length, x_p, s_prev, 'pchip');
    y_int_prev = interp1(arc_length, y_p, s_prev, 'pchip');
    if mpciter == 0
        psi_int = 0;
    else
        psi_int = atan2((y_int - y_int_prev) , (x_int - x_int_prev));
    end
    s_prev = s_0 + max((jj - 1), 0) * delta_s; 
    % scatter(mpciter, psi_int)
    % hold on

    xp0 = [xp0; x_int; y_int; psi_int; 0];
%compure ref vector
    %%used for plot
    ref = [ref; x_int, y_int, psi_int];
    
end
    xp0 = [xp0; x_int; y_int; psi_int; 0];

t0 = t0 + T;
u0 = [u(2:size(u,1),:);u(size(u,1),:)];

% up0 = reshape(u0', n_con*N,1);
up0 = [];
for jj =1:N
    up0 = [up0; u0(1,:)'];
end


% up0 = [((x_int-x_int_prev)-(y_int-y_int_prev))/(T*(cos(x0(3))+sin(x0(3)))); 
%         atan2((psi_int - psi_int_prev),(T*((x_int-x_int_prev)-(y_int-y_int_prev))...
%         /(T*(cos(x0(3))+sin(x0(3))))));
%         0];
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