clc;
clear;      % Added by Quade (just to make sure nothing is being reused)
close all;  % Uncommented by Quade
T=0.05;
t = 0:T:100;

%% Random number generator & seed 

% Below line uses a static rng seed, so results will be the same on every run
% You can use this to tune your filter by keeping outputs repeatable

% rng(0,'twister')

% Use below to "shuffle" the rng seed, so results will be different on each run
% Use this to actually test your filter

rng('shuffle')

% Below saves all information in a struct in case you need it for whatever reason

random_seed = rng;

%% FINDING TRUE DYNAMICS GIVEN F_IN %%
omega = 0.5;
f_in = 5*sin(omega*t);

x0 = [0,0]';
sol = ode45(@f_ODE_bing,t,x0);
x_True = deval(sol,t);

cd = ones(1,length(t))*0.3;
fc = ones(1,length(t));

x_True = [x_True;cd;fc];

%% DEFINING NON-LINEAR EQNS %%%
n = 4; 
m = 2;

%% PARAMETERS %%
%cd = 0.3;
%fc = 1; %update 

f = @(x,u) [
    x(1)+x(2)*T;
    x(2)+((x(3)*x(2)+x(4)*sign(x(2))))*T;
    x(3);
    x(4)
    ];

h = @(x,u)[
    x(1);
    (x(3)*x(2)+x(4)*sign(x(2)))
];

x = zeros(n, length(t)); % Initializes states to zero
Q = 1e-6*eye(n); 
R = 1e-3*eye(m); 

x(:,1) = [0,0,0,0];
x_ukf = x;
P_ukf = 10*Q;
w = mvnrnd(zeros(length(t),n),Q)'; 
v = mvnrnd(zeros(length(t),m),R)';

%% Added by Quade - adding process noise
x_True = x_True + w;

f_measured = zeros(m,length(t));
%% UNSCENTED KALMAN FILTER%%
for k = 1:length(t)-1
    f_measured(:,k) = h(x_True(:,k),0) + v(k);
    [x_ukf(:,k+1), P_ukf(:,:,k+1)] = ukf_input(x_ukf(:,k), f_measured(:,k),0, P_ukf(:,:,k), f, h, Q, R); % Calls UKF function
end

SumE = zeros(n,1);

rmse = sqrt( mean(x_True - x_ukf,2).^2 ) % this is just vectorized
% rmse_force = sqrt(mean(f_measured - f_in,2).^2 ) % this is just vectorized

tiledlayout(2,2);
nexttile;
plot(t, x_ukf(1,:)); hold all; plot(t,x_True(1,:));title("Position"); legend("UKF","True");
xlabel('Time(s)')
ylabel('x')
nexttile;
plot(t, x_ukf(2,:)); hold all; plot(t,x_True(2,:));title("Velocity"); legend("UKF","True");
xlabel('Time(s)')
ylabel('v')
nexttile;
plot(t,x_ukf(3,:)); hold all; plot(t,x_True(3,:));title("Cd"); legend("UKF","True");
ylim([-2,3])
xlabel('Time(s)')
ylabel('cd')
nexttile;
plot(t,x_ukf(4,:)); hold all; plot(t,x_True(4,:));title("Fc"); legend("UKF","True");
xlabel('Time(s)')
ylabel('fc')
ylim([-2,3])

figure;
plot(t, f_measured); 
hold on 
% plot(x_True(2,:),f_measured);
title("Force Velocity Hysteresis");
legend("Measured","Input");
xlabel('V')
ylabel('F')
hold off

% figure;plot(x_ukf(2,:),f_in-f_measured); hold all; plot(x_True(2,:),f_in-f_measured);title("FV"); legend("UKF","True");
% figure; plot(t, x_ukf(2,:)); hold all; plot(t,x_True(2,:));title("Velocity"); legend("UKF","True");
% figure; plot(t,f_measured); hold all; plot(t, f_in); title("Force"); legend("Damper force", "Actuator Force");

