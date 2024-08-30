clc;clear;
% close all;

%% PARAMETERS %%
T=0.05;
t = 0:T:30;

n = 3; 
m = 1;

omega = 1;
cd = 1200;
fc = 30;
f0 = -7; 

%% TRUE DYNAMICS %%
period = 15;
omega = 2*pi*1/period;
x = -0.0375*sin(omega*t-pi/2)-0.0375;

xdot = [diff(x)/T 0];
% xdot = 5*omega*cos(omega*t);
% xdot = sawtooth(omega * t, 0.5);;

cd = ones(1,length(t))*cd;
fc = ones(1,length(t))*fc;
f0 = ones(1,length(t))*f0;

x_True = [cd' fc' f0']';

% STATE VECTOR
f = @(x,u)[
    x(1)
    x(2)
    x(3)
    ];

% MEASUREMENT VECTOR
h = @(x,u)[
    x(1)*100*u+x(2)*sign(u)+x(3);
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
    f_measured(k) = h(x_True(:,k),xdot(k))+v(k);
    [x_ukf(:,k+1), P_ukf(:,:,k+1)] = ukf_input(x_ukf(:,k), f_measured(k),xdot(k), P_ukf(:,:,k), f, h, Q, R); % Calls UKF function
end

steady = 25;
indic = find(t == 25);
f_ukf = zeros(m,length(t));

for k = 1:length(t)-1
    f_ukf(k) = x_ukf(1,indic)*100*xdot(k)+x_ukf(2,indic)*sign(xdot(k))+x_ukf(3,indic);
end

tiledlayout(2,2);
nexttile;
plot(t, x_ukf(1,:)); hold all; plot(t,cd); ylim([0,1500]);title("c0 Parameter ID"); legend("UKF","True"); xlabel("Time(s)"); ylabel("Ns/m");
nexttile;
plot(t, x_ukf(2,:)); hold all; plot(t,fc);title("fc Parameter ID"); legend("UKF","True");xlabel("Time(s)"); ylabel("N");
nexttile;
plot(t, x_ukf(3,:)); hold all; plot(t,f0);title("f0 Parameter ID"); legend("UKF","True");xlabel("Time(s)"); ylabel("N");
nexttile;
plot(t, f_measured); hold all; plot(t,f_ukf);title("Simulated Measured vs. UKF"); legend("UKF","True");xlabel("Time(s)"); ylabel("N");


