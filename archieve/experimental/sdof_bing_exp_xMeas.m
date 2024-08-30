clc;
close all;
clear
T = 0.0147;

tic;
data = readtable(['data_voltagevarying_0volt2.csv']);

t = data(:,2);
t = t{:,:};
f_measured = data(:,6);
f_measured= f_measured{:,:};

x0 = data(:,3);
x0 = x0{:,:};
x1 = data(:,4);
x1 = x1{:,:};

% Filter velocity
order = 3;
framelen = 11;
x1 = sgolayfilt(x1,order,framelen);

% x1 = diff(x0)/T;
% x1 = [diff(x0)/T;0];
% size(x1)
% size(x0)
x_True = [x0 x1]';

%% STATE AND MEASUREMENT VECTOR SIZE %%%
n = 5; 
m = 1;

%% PARAMETERS from LLS %%
cd= 271.56716265121725;
fc= 5.529491411745272;
f0= -8.055519974161356;
% cd = 352.5194;   
% fc = 5.5396;
% f0 = -8.2707;

%% STATE VECTOR %%
f = @(x,u) [
    x(1); %fc estimate
    x(2); %f0 estimate
    x(3); %cd estimate
    ];

h = @(x,u)[
    u(1)
    u(2)
    x(5)*u(2)+x(3)*sign(u(2))+x(4);
];

x_ukf = zeros(n, length(t)); % Initializes states to zero

x_ukf(:,1) = [0,0,0,0,0]; %inital values

Q = 1e-2*eye(n); 
R = 1e-2*eye(m); 
P_ukf = 10*Q;

%% UNSCENTED KALMAN FILTER%%
for k = 1:length(t)-1
    [x_ukf(:,k+1), P_ukf(:,:,k+1)] = ukf_input(x_ukf(:,k), f_measured(k), x_True(:,k), P_ukf(:,:,k), f, h, Q, R); % Calls UKF function
end

%% FIND AVG. OF PARAM FROM UKF %%
fc_mean = mean(x_ukf(3,:));
f0_mean = mean(x_ukf(4,:));
cd_mean = mean(x_ukf(5,:));

f_mean = zeros(1, length(t)); % Initializes measurements to zero
f_lls = zeros(1, length(t)); % Initializes measurements to zero

for k = 1:length(t)-1
    f_mean(k) = cd_mean*x_True(2,k)+fc_mean*sign(x_True(2,k))+f0_mean;
    f_lls(k) = cd*x_True(2,k)+fc*sign(x_True(2,k))+f0;
end

elapsed = toc

squaredDifferences = (f_measured' - f_lls).^2;
mse = mean(squaredDifferences);
rmse_lls = sqrt(mse);

squaredDifferences = (f_measured' - f_mean).^2;
mse = mean(squaredDifferences);
rmse_mean = sqrt(mse);

paramName = {'Fc', 'F0', 'Cd'};
paramUKF = [fc_mean, f0_mean, cd_mean];
paramLLS = [fc, f0, cd];
Tparam = table(paramName', paramUKF', paramLLS', 'VariableNames', {'Parameter', 'UKF', 'LLS'});
disp(Tparam);

method = {'UKF', 'LLS'};
rmseCol = [rmse_mean, rmse_lls];
Trmse = table(method', rmseCol', 'VariableNames', {'Method', 'RMSE'});
disp(Trmse);

tiledlayout(2,3);
nexttile;
plot(t, x_ukf(1,:)); hold all; plot(t,x_True(1,:));title("Position"); legend("UKF","True");
xlabel('Time(s)')
ylabel('x')
nexttile;
plot(t, x_ukf(2,:)); hold all; plot(t,x_True(2,:));title("Velocity"); legend("UKF","True");
xlabel('Time(s)')
ylabel('v')
nexttile;
plot(t, x_ukf(3,:)); hold all; title("Fc"); legend("UKF");
xlabel('Time(s)')
ylabel('Fc')
nexttile;
plot(t, x_ukf(4,:)); hold all; title("F0"); legend("UKF");
xlabel('Time(s)')
ylabel('F0')
nexttile;
plot(t, x_ukf(5,:)); hold all; title("Cd"); legend("UKF");
xlabel('Time(s)')
ylabel('Cd')
nexttile;
plot(t, f_measured); 
hold all;
plot(t,f_mean);
plot(t,f_lls);
title("Force Vs. Time Post Processing");
legend("Measured","UKF","LLS");
xlabel('Time(s)')
ylabel('F')

