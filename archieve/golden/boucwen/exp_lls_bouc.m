clc;
close all;
clear;

tic;

% Load or generate your dataset
data = readtable(['v250000_a50_121mv2_clean.csv']);

t = data(:,2);
t = t{:,:};
f_measured = data(:,6);
f_measured= f_measured{:,:};
x1 = data(:,4);
x1 = x1{:,:};

order = 3;
framelen = 11;
x1 = sgolayfilt(x1,order,framelen)*100;

x2 = data(:,3);
x2 = x2{:,:}*100;

c0= 8.86846727390849;
k0= 0.597632140788303;
alpha= 30.837990991911244;
beta= -954.4497438208226;
gamma= 972.7733104878047;
A= 7.622000649825172;
x0= 5.388211694664788;
n_po = 2;

z = zeros(1,length(t)-1)';
z(1) = 0.1;
T = 0.0147;
for k = 1:length(t)-1
    z(k+1) = z(k)+(-gamma*abs(x1(k))*z(k)*abs(z(k))^(n_po-1)-beta*x1(k)*abs(z(k))^n_po+A*x1(k))*T;
end

% Define the parametric model
bw = @(p, vel, pos, evo) p(1)*vel+p(2)*(pos-p(3))+p(4)*evo;

% Initial guess for the coefficients [amplitude, center, width]
initial_guess = [0, 0, 0, 0];

% Define the objective function (residuals)
obj_function = @(p) f_measured - bw(p,x1,x2,z);

% Perform nonlinear least squares optimization
options = optimset('Display', 'iter'); % Display optimization process
coefficients = lsqnonlin(obj_function, initial_guess, [], [], options);

% Generate fitted curve using the optimized coefficients
fitted_curve = bw(coefficients, x1, x2,z);

toc
size(f_measured)
size(fitted_curve)
squaredDifferences = (f_measured - fitted_curve).^2;
% Compute the mean of squared differences
mse = mean(squaredDifferences);
% Compute the RMSE by taking the square root of MSE
rmse = sqrt(mse)

% rmse_force = sqrt(mean(f_measured' - fitted_curve',2).^2)% this is just vectorized


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
