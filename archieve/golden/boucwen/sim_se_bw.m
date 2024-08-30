clc;clear;close all;

T=0.05;
t = 0:T:100;

n = 3; 
m = 1;

%% TRUE DYNAMICS %%
x0 = [0,0,0]';
sol = ode45(@boucwen,t,x0);
x_ODE = deval(sol,t);
z = x_ODE(3,:);

period = 15;
w = 2*pi*1/period;
x = (-0.0375*sin(w*t-pi/2)-0.0375)*100;
xdot = -0.0375*w*cos(w*t-pi/2)*100;
xddot = 0.0375*w^2*sin(w*t-pi/2)*100;
x_True = [x' xdot' z']' ;

%% PARAMETERS %%
c0= 11.825884573691747
k0= 0.6168415593308919
alpha= 44.45425042190868
beta= -3061.832988105183
gamma= 3078.384741848536
A= 2.761665766402648
x0= 7.191401303755131
n_po = 2;

% STATE VECTOR
f = @(x,u)[
    x(1)+x(2)*T;
    x(2)+(u-(c0*x(2)+k0*(x(1)-x0)+alpha*x(3)))*T;
    x(3)+(-gamma*abs(x(2))*x(3)*abs(x(3))^(n_po-1)-beta*x(2)*abs(x(3))^n_po+A*x(2))*T;
];

% MEASUREMENT VECTOR
h = @(x,u)[
    c0*x(2)+k0*(x(1)-x0)+alpha*x(3)
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
f_ukf = zeros(m,length(t));

%% UNSCENTED KALMAN FILTER%%
for k = 1:length(t)-1
    f_measured(k) = -h(x_True(:,k),0)+v(k);
    f_ukf(k) = c0*x_ukf(2,k)+k0*(x_ukf(1,k)-x0)+alpha*x_ukf(3,k);
    [x_ukf(:,k+1), P_ukf(:,:,k+1)] = ukf_input(x_ukf(:,k), f_measured(k),xddot(k), P_ukf(:,:,k), f, h, Q, R); % Calls UKF function
end

tiledlayout(2,2)
nexttile;
plot(t,x_ukf(1,:))
hold on
plot(t,x_True(1,:))
legend("UKF", "True")
title("Position estimation")
xlabel("Time(s)"); ylabel("x (cm)");
nexttile;
plot(t,x_ukf(2,:))
hold on
plot(t,x_True(2,:))
legend("UKF", "True")
title("Velocity estimation")
xlabel("Time(s)"); ylabel("v(cm/s)");
nexttile;
plot(t,x_ukf(3,:))
hold on
plot(t,x_True(3,:))
legend("UKF", "True")
title("Evloutionary Variable")
xlabel("Time(s)"); ylabel("");
nexttile;
plot(t,f_ukf)
hold on 
plot(t,f_measured)
legend("UKF", "True")
title("Evloutionary Variable")
xlabel("Time(s)"); ylabel("N");
