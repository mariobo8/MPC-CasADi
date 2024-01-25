clearvars;
close all
clc


addpath('~/casadi-3.6.4-linux64-matlab2018b')
line_width = 1.5;
fontsize_labels = 14;
import casadi.*
theta = SX.sym('theta');
mu = [6; 20; 5; 0.35];  
rho = 15*sin(theta * mu(4)) .* (1/(1+abs(theta)));
f2 = Function('f',{theta},{rho});
theta_range = linspace(-30, 0, 300);
rho_val = full(f2(theta_range)); 
%plot path
    plot(theta_range, rho_val,'linewidth',1.5);
        ylabel('$\rho(\theta)$ (m)','interpreter','latex','FontSize',fontsize_labels)
    xlabel('$\theta$ (m)','interpreter','latex','FontSize',fontsize_labels)

    axis([-30 0 -2.5 5.5]) 
        box on;
    grid on