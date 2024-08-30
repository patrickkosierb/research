clc;
clear;
close all;
T=0.05;
t = 0:T:100;

factor = 2;

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
% f_in = ones(1,length(t))*5;
x0 = [0,0,0]';
sol = ode45(@f_ODE,t,x0);
x_True = deval(sol,t);

c = ones(1,length(t))*0.3*factor;
k = ones(1,length(t))*1;

x_True = [x_True;c;k];

%% DEFINING NON-LINEAR EQNS %%%
n = 5; 
m = 1;

%c = 0.6;
%k = 1;
beta = 4;
gamma = 1;
n_pa = 2;
% k=1;

f = @(x,u) [
    x(1)+x(2)*T;
    x(2)+(u-(x(4)*x(2)+x(5)*x(3)))*T;
    x(3)+(x(2)-gamma*abs(x(2))*x(3)*abs(x(3))^(n_pa-1)-beta*x(2)*abs(x(3))^n_pa)*T;
    x(4);
    x(5)
    ];

h = @(x,u)[
    (x(4)*x(2)+x(5)*x(3)) 
];

x = zeros(n, length(t)); % Initializes states to zero
Q = 1e-6*eye(n); 
R = 1e-3*eye(m); 
x(:,1) = [0,0,0,0,0]';
x_ukf = x;
P_ukf = 10*Q;
w = mvnrnd(zeros(length(t),n),Q)'; 
v = mvnrnd(zeros(length(t),m),R)'; 

%% Added by Quade - adding process noise
x_True = x_True + w;

%% CREATING RANDOM F_MR DATA %%
f_measured = zeros(1,length(t));
f_ukf = zeros(1,length(t));

%% UNSCENTED KALMAN FILTER%%
for k = 1:length(t)-1
    f_measured(k) = h(x_True(:,k),f_in(k))+v(k);
    f_ukf(k) = h(x_ukf(:,k),f_in(k));
    [x_ukf(:,k+1), P_ukf(:,:,k+1)] = ukf_input(x_ukf(:,k), f_measured(k),f_in(k), P_ukf(:,:,k), f, h, Q, R); % Calls UKF function
end

SumE = zeros(n,1);
Fdiff= zeros(1,1);

rmse_state = sqrt(mean(x_True - x_ukf,2).^2 ) % this is just vectorized
% rmse_force = sqrt(mean(f_measured - f_in,2).^2 ) % this is just vectorized


% figure; plot(t, x_ukf(1,:)); hold all; plot(t,x_True(1,:));title("Position"); legend("UKF","True");
% figure; plot(t, x_ukf(2,:)); hold all; plot(t,x_True(2,:));title("Velocity"); legend("UKF","True");
% figure; plot(t,f_measured); hold all; plot(t, f_in); title("Force"); legend("Damper force", "Actuator Force");
% figure;
tiledlayout(2,2);
nexttile
plot(t, x_ukf(1,:)); 
hold on 
plot(t,x_True(1,:));
title("Position"); 
legend("UKF","True");
xlabel('Time(s)')
ylabel('x')
% ylim ([-1 5])
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
plot(t, x_ukf(4,:)); 
hold on 
plot(t,x_True(4,:));
title("C");
legend("UKF","True");
xlabel('Time(s)')
ylim ([-1 1])
hold off

nexttile
plot(t, x_ukf(5,:)); 
hold on 
plot(t,x_True(5,:));
title("K");
legend("UKF","True");
ylim ([-2 2*factor])
xlabel('Time(s)')
hold off

% nexttile
figure;
plot(x_ukf(2,:), f_in-f_measured); 
hold on 
plot(x_True(2,:),f_in-f_measured);
title("Force Velocity Hysteresis");
legend("Measured","Input");
xlabel('V')
ylabel('F')
hold off

% nexttile
figure;
plot(t, f_measured); 
hold on 
plot(t,f_ukf);
title("Force vs. Time");
legend("True","UKF");
xlabel('t')
ylabel('F')
hold off
