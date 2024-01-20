clearvars;
close all;
clc;
theta_range = linspace(-30, 0, 300);  % Create a time vector from 0 to 2*pi

% Parametric equations for x and y coordinates
alpha = 6; gamma = 20; beta = 5; w = 0.35;  
rho = @(theta) - alpha * log(gamma./(beta + abs(theta))) .* sin(w * theta);  % Adjust the amplitude as needed


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


