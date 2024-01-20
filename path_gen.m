clearvars;
close all;
clc;
theta_range = linspace(0, 2*pi, 300);  % Create a time vector from 0 to 2*pi

% Parametric equations for x and y coordinates

rho = @(theta) 1./(1+5*exp(-theta)).*sin(theta) + 10;  % Adjust the amplitude as needed


% psi = diff(y)./diff(x);
% psi = [psi(1), psi];

rho_val = rho(theta_range);

% Plot the function
figure;
plot(theta_range, rho_val, 'LineWidth', 2);
title('\rho(\theta) Plot');
xlabel('\theta');
ylabel('\rho(\theta)');
grid on;


