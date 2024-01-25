clc; clearvars; close all;

load("errore.mat");
load("average_time.mat");

for jj = 1 : 7 
    for kk = 1 : 7
        for nn = 1 : 7
            for ll = 1 : 4
                e_m.x(jj).y_m(kk).psi_m(nn).theta_m(ll,:) = ...
                mean(e.x(jj,:).y(kk,:).psi(nn,:).theta(ll,:));
            end
        end
    end
end
min_val = 10;
for jj = 1 : 7 
    for kk = 1 : 7
        for nn = 1 : 7
            for ll = 1 : 4
                if e_m.x(jj).y_m(kk).psi_m(nn).theta_m(ll) < min_val
                    a = jj; b = kk; c = nn; d = ll;
                    min_val = e_m.x(jj).y_m(kk).psi_m(nn).theta_m(ll);
                end
            end
        end
    end
end
e_val = e_m.x(a).y_m(b).psi_m(c).theta_m(d) 