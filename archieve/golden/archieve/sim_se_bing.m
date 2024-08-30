clc;clear;close all;

%% PARAMETERS %%
T=0.05;
t = 0:T:30;

n = 2; 
m = 1;

omega = 1;
cd = 1200;
fc = 30;
f0 = -7; 

%% TRUE DYNAMICS %%
x = 10*sin(omega*t);
xdot = 10*omega*cos(omega*t);
xddot = -10*omega^2*sin(omega*t);
x_True = [x;xdot];

% STATE VECTOR
f = @(x,u)[
    x(1)+x(2)*T;
    x(2)+u*T;
    ];

% MEASUREMENT VECTOR
h = @(x,u)[
    cd*x(2)+fc*sign(x(2));
];

% NOISE
Q = 1e-2*eye(n); 
R = 1e-2*eye(m); 
P_ukf = 10*Q;

rng('shuffle')
random_seed = rng;
w = mvnrnd(zeros(length(t),n),Q)'; 
v = mvnrnd(zeros(length(t),m),R)'; 

x_True = x_True + w;
f_measured = zeros(m,length(t));
x_ukf = zeros(n,length(t));

%% UNSCENTED KALMAN FILTER%%
for k = 1:length(t)-1
    f_measured(k) = h(x_True(:,k),xddot(k))+v(k);
    [x_ukf(:,k+1), P_ukf(:,:,k+1)] = ukf_input(x_ukf(:,k), f_measured(k),xddot(k), P_ukf(:,:,k), f, h, Q, R); % Calls UKF function
end

% rmse = sqrt( mean(x_True - x_ukf,2).^2 ) % this is just vectorized
% rmse_force = sqrt(mean(f_measured - f_in,2).^2 ) % this is just vectorized


tiledlayout(1,2);
nexttile;
plot(t, x_ukf(1,:)); hold all; plot(t,x_True(1,:));title("Position"); legend("UKF","True"); xlabel("Time(s)"); ylabel("x");
nexttile;
plot(t, x_ukf(2,:)); hold all; plot(t,x_True(2,:));title("Velocity"); legend("UKF","True");xlabel("Time(s)"); ylabel("v");
% % nexttile;
% % plot(x_ukf(2,:),f_in-f_measured); hold all; plot(x_True(2,:),f_in-f_measured);title("FV"); legend("UKF","True");


% % figure; plot(t, x_ukf(2,:)); hold all; plot(t,x_True(2,:));title("Velocity"); legend("UKF","True");
% figure; plot(t,f_measured); hold all; plot(t, f_in); title("Force"); legend("Damper force", "Actuator Force");

