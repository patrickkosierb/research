clc;
clear;
close all;

T = 0.0147;
data = readtable(['data_voltagevarying_0volt.csv']);

t = data(:,2);
t = t{:,:};
f_measured = data(:,6);
f_measured= f_measured{:,:};

x0 = data(:,3);
x0 = x0{:,:};
x1 = data(:,4);
x1 = x1{:,:};

order = 3;
framelen = 11;
x1 = sgolayfilt(x1,order,framelen);

x_True = [x0 x1]';

omega = 0.15;
f_in = 10*sin(omega*t)-7;

%% DEFINING NON-LINEAR EQNS %%%
n = 3; 
m = 1;

gamma = 130.79645522568575;
beta = -0.89773144490448;
n_pa = 2;
k_var = -60.078547581405964;
c= 1201.3445613819722;
mass=0.5

f = @(x,u) [
    x(1)+x(2)*T;
    x(2)+(c*x(2)+k_var*x(3)-u)*T/mass;
    x(3)+(x(2)-beta*abs(x(2))*x(3)*abs(x(3))^(n_pa-1)-gamma*x(2)*abs(x(3))^n_pa)*T;
    ];

h = @(x,u)[
    c*x(2)+k_var*x(3);
];

Q = 1e-2*eye(n); 
R = 1e-2*eye(m); 

x_ukf = zeros(n, length(t));
f_ukf = zeros(1,length(t));

x_ukf(:,1) = [0,0,0.1]';

P_ukf = 10*Q;

%% UNSCENTED KALMAN FILTER%%
for k = 1:length(t)-1
    f_ukf(k) = h(x_ukf(:,k),x_True(1:2,k));
    [x_ukf(:,k+1), P_ukf(:,:,k+1)] = ukf_input(x_ukf(:,k), f_measured(k),f_in(k), P_ukf(:,:,k), f, h, Q, R); % Calls UKF function
end


% squaredDifferences = (x_True(1:2,:) - x_ukf(1:2,:)).^2;
% mse = mean(squaredDifferences);
% rmse_mean = sqrt(mse);


tiledlayout(2,2);
nexttile
plot(t, x_ukf(1,:)); 
hold on 
plot(t,x_True(1,:));
title("Position"); 
legend("UKF","True");
xlabel('Time(s)')
ylabel('x')
hold off
nexttile
plot(t, x_ukf(2,:)); 
hold on 
plot(t,x_True(2,:));
title("Velocity"); 
legend("UKF","True");
xlabel('Time(s)')
ylabel('v')
hold off
nexttile
plot(t, x_ukf(3,:)); 
hold on 
title("Z");
legend("UKF","True");
xlabel('Time(s)')
hold off

figure;
plot(t, f_measured); 
hold on 
plot(t, f_in); 
title("Force vs. Time");
legend("Measured","Mean");
xlabel('t')
ylabel('F')
hold off

