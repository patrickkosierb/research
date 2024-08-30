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

x0 = [0,0,0]';
sol = ode45(@f_ODE,t,x0);
x_True = deval(sol,t);

factor = 2;
c = ones(1,length(t))*0.3*factor;
k = ones(1,length(t))*1;
beta = ones(1,length(t))*2*factor;
gamma = ones(1,length(t))*1;

x_True = [x_True;c;k;beta;gamma];

%% DEFINING NON-LINEAR EQNS %%%
n = 7; 
m = 1;

% c=0.3;
% k = 1;
% beta = 2;
% gamma = 1;
n_pa = 2;

f = @(x,u) [
    x(1)+x(2)*T;
    x(2)+(u-(x(4)*x(2)+x(5)*x(3)))*T;
    x(3)+(x(2)-x(6)*abs(x(2))*power(abs(x(3)),n_pa-1)*x(3)-x(7)*x(2)*abs(x(3))^n_pa)*T;
    x(4);
    x(5);
    x(6);
    x(7);
    ];

h = @(x,u)[
    u-(x(4)*x(2)+x(5)*x(3))
];

x = zeros(n, length(t)); % Initializes states to zero
Q = 1e-6*eye(n); 
R = 1e-3*eye(m); 

% x(:,1) = [0,0,0,0.2,0.7,1.5]';
x(:,1) = [0,0,0,0.1*factor,0.1,1*factor,0.1*factor]';
x_ukf = x;
P_ukf = 10*Q;
w = mvnrnd(zeros(length(t),n),Q)'; 
v = mvnrnd(zeros(length(t),m),R)'; 

%% Added by Quade - adding process noise
x_True = x_True + w;

%% CREATING RANDOM F_MR DATA %%
f_measured = zeros(1,length(t));

%% UNSCENTED KALMAN FILTER%%
for k = 1:length(t)-1
    f_measured(k) = h(x_True(:,k),f_in(k))+v(k);
    [x_ukf(:,k+1), P_ukf(:,:,k+1)] = ukf_input2(x_ukf(:,k), f_measured(k),f_in(k), P_ukf(:,:,k), f, h, Q, R); % Calls UKF function
end

SumE = zeros(n,1);
Fdiff= zeros(1,1);

rmse_state = sqrt(mean(x_True - x_ukf,2).^2 ) % this is just vectorized
rmse_force = sqrt(mean(f_measured - f_in,2).^2 ) % this is just vectorized

% figure;
tiledlayout(2,3);
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

% nexttile
% plot(t, x_ukf(3,:)); 
% hold on 
% plot(t,x_True(3,:));
% title("Z"); 
% legend("UKF","True");
% hold off

% nexttile
% plot(x_ukf(1,:),x_ukf(3,:));
% hold on;
% plot(x_True(1,:),x_True(3,:))
% title("X Vs. Z");
% legend("UKF","True");
% hold off;

% nexttile
% plot(t, f_measured); 
% hold on 
% plot(t,f_in);
% title("F");
% legend("Measured","Input");
% hold off

nexttile
plot(t, x_ukf(4,:)); 
hold on 
plot(t,x_True(4,:));
title("C");
legend("UKF","True");
xlabel('Time(s)')
% ylabel('')
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

nexttile
plot(t, x_ukf(6,:)); 
hold on 
plot(t,x_True(6,:));
title("Beta");
legend("UKF","True");
ylim ([-1 3*factor])
xlabel('Time(s)')
hold off

nexttile
plot(t, x_ukf(7,:)); 
hold on 
plot(t,x_True(7,:));
title("Gamma");
legend("UKF","True");
ylim ([-1 2*factor])
xlabel('Time(s)')
hold off
% 
% figure;
% plot(x_ukf(1,:), x_ukf(3,:)); 
% hold on 
% plot(x_True(1,:), x_True(3,:)); 
% title("Z vs. X");
% % ylim ([-0.8 0.8])
% legend("UKF","True");
% hold off

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




