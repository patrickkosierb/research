clc;clear;close all;

T=0.05;
t = 0:T:100;

n = 1; 
m = 1;

%% TRUE DYNAMICS %%
x0 = [0,0,0]';
sol = ode45(@boucwen,t,x0);
x_ODE = deval(sol,t);
z = x_ODE(3,:);
x_True = [z']';

%% INPUT %%
period = 15;
w = 2*pi*1/period;
x = (-0.0375*sin(w*t-pi/2)-0.0375)*100;
xdot = -0.0375*w*cos(w*t-pi/2)*100;
xddot = 0.0375*w^2*sin(w*t-pi/2)*100;

in = [x' xdot']';

%% PARAMETERS %%
c0= 11.825884573691747;
k0= 0.6168415593308919;
alpha= 44.45425042190868;
beta= -3061.832988105183;
gamma= 3078.384741848536;
A= 2.761665766402648;
x0= 7.191401303755131;
n_po = 2;

% STATE VECTOR
f = @(z,u)[
    z+(-gamma*abs(u(2))*z*abs(z)^(n_po-1)-beta*u(2)*abs(z)^n_po+A*u(2))*T;
];

% MEASUREMENT VECTOR
h = @(z,u)[
    c0*u(2)+k0*(u(1)-x0)+alpha*z
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

for k = 1:length(t)-1
    f_measured(k) = h(x_True(:,k),in(:,k))+v(k);
    [x_ukf(:,k+1), P_ukf(:,:,k+1)] = ukf_input(x_ukf(:,k), f_measured(k),in(:,k), P_ukf(:,:,k), f, h, Q, R); % Calls UKF function
end


plot(t,x_ukf(1,:))
hold on
plot(t,x_True(1,:))
