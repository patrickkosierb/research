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


%% DEFINING NON-LINEAR EQNS %%%
n = 5; 
m = 1;

beta = 130.79645522568575;
gamma = -0.89773144490448;
n_pa = 2;
k_var = -60.078547581405964;
c= 1201.3445613819722;
mass=0.5

f = @(x,u) [
    u(1);
    u(2);
    x(3)+(u(2)-gamma*abs(u(2))*x(3)*abs(x(3))^(n_pa-1)-x(4)*u(2)*abs(x(3))^n_pa)*T;
    x(4);
    x(5)
    ];

h = @(x,u)[
    x(5)*u(2)+k_var*x(3);
];

Q = 1e-2*eye(n); 
R = 1e-2*eye(m); 

x_ukf = zeros(n, length(t));
f_ukf = zeros(1,length(t));

x_ukf(:,1) = [0,0,0.1,0,0]';

P_ukf = 10*Q;

%% UNSCENTED KALMAN FILTER%%
for k = 1:length(t)-1
    [x_ukf(:,k+1), P_ukf(:,:,k+1)] = ukf_input(x_ukf(:,k), f_measured(k),x_True(1:2,k), P_ukf(:,:,k), f, h, Q, R); % Calls UKF function
end

bmean =  mean(x_ukf(4,:));
cmean =  mean(x_ukf(5,:));

f_mean = zeros(1,length(t));
z_mean = zeros(1,length(t));

f_true = zeros(1,length(t));
z_true = zeros(1,length(t));

for k = 1:length(t)-1
    z_true(k+1) = z_true(k)+(x1(k)-gamma*abs(x1(k))*z_true(k)*abs(z_true(k))^(n_pa-1)-beta*x1(k)*abs(z_true(k))^n_pa)*T;
    f_true(k+1) = c*x1(k)+k_var*z_true(k+1);

    z_mean(k+1) = z_mean(k)+(x1(k)-gamma*abs(x1(k))*z_mean(k)*abs(z_mean(k))^(n_pa-1)-bmean*x1(k)*abs(z_mean(k))^n_pa)*T;
    f_mean(k+1) = cmean*x1(k)+k_var*z_mean(k);
end

squaredDifferences = (f_measured' - f_mean).^2;
mse = mean(squaredDifferences);
rmse_mean = sqrt(mse);


tiledlayout(2,3);
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
nexttile
plot(t, x_ukf(4,:)); 
hold on 
title("Beta");
legend("UKF","True");
xlabel('Time(s)')
hold off
nexttile
plot(t, x_ukf(5,:)); 
hold on 
title("C");
legend("UKF","True");
xlabel('Time(s)')
hold off
nexttile
plot(t, f_measured); 
hold on 
plot(t,f_mean);
hold on 
plot(t,f_true);
title("Force vs. Time");
legend("Measured","UKF","LLS");
xlabel('t')
ylabel('F')
hold off
