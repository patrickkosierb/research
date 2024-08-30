clc;
close all;
clear;

tic;

% Load or generate your dataset
data = readtable(['data_voltagevarying_0volt.csv']);

t = data(:,2);
t = t{:,:};
f_measured = data(:,6);
f_measured= f_measured{:,:};
x1 = data(:,4);
x1 = x1{:,:};

order = 3;
framelen = 11;
x1 = sgolayfilt(x1,order,framelen);


% Define the parametric model
% x(3)+*T;
n_pa = 2

bouc = @(p, x) p(1)*x+p(2)*integral((x(2)-p(3)*abs(x(2))*x(3)*abs(x(3))^(n_pa-1)-p(4)*x(2)*abs(x(3))^n_pa),0,inf);

% Initial guess for the coefficients [amplitude, center, width]
initial_guess = [0, 0, 0];

% Define the objective function (residuals)
obj_function = @(p) f_measured - bouc(p, x1);

% Perform nonlinear least squares optimization
options = optimset('Display', 'iter'); % Display optimization process
coefficients = lsqnonlin(obj_function, initial_guess, [], [], options);

% Generate fitted curve using the optimized coefficients
fitted_curve = bouc(coefficients, x1);

squaredDifferences = (f_measured - fitted_curve).^2;
% Compute the mean of squared differences
mse = mean(squaredDifferences);
% Compute the RMSE by taking the square root of MSE
rmse = sqrt(mse)

% rmse_force = sqrt(mean(f_measured' - fitted_curve',2).^2)% this is just vectorized


elapsed = toc
% Plot the results
figure;
plot(t, f_measured, 'b', 'MarkerSize', 3); % Plot data points
hold on;
plot(t, fitted_curve, 'r-', 'LineWidth', 2); % Plot fitted curve
xlabel('x');
ylabel('y');
title('Nonlinear Least Squares Fitting');
legend('Data', 'Fitted Curve');

% Display the optimized coefficients
disp('Optimized coefficients:');
disp(coefficients);
