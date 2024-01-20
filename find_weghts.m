clearvars;
close all
clc
e = struct;
avg_time = struct;
for jj = 1 : 7 
    for kk = 1 : 7
        for nn = 1 : 7
            for ll = 1 : 4
                Q = zeros(4,4); Q(1,1) = 1*10^(jj); Q(2,2) = 1*10^(kk); 
                                Q(3,3) = 1*10^(nn); q(4,4) = 1*10^(ll);% weighing matrices (states)
                R = zeros(3,3); R(1,1) = 10; R(2,2) = 10; ...
                                R(3,3) = 1;% weighing matrices (controls)
                eps = 1500;
                run('carexample_path_tracking_MS.m')
                e.x(jj,:).y(kk,:).psi(nn,:).theta(ll,:) = error;
                avg_time.x(jj).y(kk).psi(nn).theta(ll) = average_mpc_time;
            end
        end
    end
end