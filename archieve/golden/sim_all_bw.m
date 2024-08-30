clc;clear;close all;

T=0.05;
t = 0:T:30;

n = 8; 
m = 1;

omega = 1;

%% TRUE DYNAMICS %%
x = 5*sin(omega*t);
xdot = 5*omega*cos(omega*t);
in = [x' xdot']';

x0 = [0,0,0]';
sol = ode45(@boucwen,t,x0);
x_ODE = deval(sol,t);
z =x_ODE(3,:);

%% PARAMETERS %%
c0= 7.0275110546299855;
k0= 0.02110224311505027;
alpha= 0.012528469307446598;
beta= -0.6284883734431349;
gamma= 0.6277882786656479;
A= 39.30100479335722;
x0= 379.45901441054866;
n_po = 2;

c0= ones(1,length(t))*c0;
k0= ones(1,length(t))*k0;
alpha= ones(1,length(t))*alpha;
beta= ones(1,length(t))*beta;
gamma= ones(1,length(t))*gamma;
A= ones(1,length(t))*A;
x0= ones(1,length(t))*x0;
n_po = 2;

x_True = [z' c0' k0' alpha' beta' gamma' A' x0']';

% STATE VECTOR
f = @(x,u)[
    x(1) + (-x(6)*abs(u(2))*x(1)*abs(x(1))^(n_po-1)-x(5)*u(2)*abs(x(1))^n_po+x(7)*u(2))*T;
    x(2); %c0
    x(3); %k0
    x(4); %alpha
    x(5); %beta
    x(6); %gamma
    x(7); %A
    x(8); %x0
    ];

% MEASUREMENT VECTOR
h = @(x,u)[
    x(2)*u(2)+x(3)*(u(1)-x(8))+x(4)*x(1);
];

% NOISE
Q = 1e-2*eye(n); 
R = 1e-2*eye(m); 
P_ukf = 10*Q;

rng('shuffle')
random_seed = rng;
w = mvnrnd(zeros(length(t),n),Q)'; 
v = mvnrnd(zeros(length(t),m),R)'; 

% x_True = x_True + w;
f_measured = zeros(m,length(t));
x_ukf = zeros(n,length(t));
x_ukf(:,1) = [0.5,1,1,1,1,1,1,1]';

%% UNSCENTED KALMAN FILTER%%
for k = 1:length(t)-1
    f_measured(k) = h(x_True(:,k),in(:,k))+v(k);
    [x_ukf(:,k+1), P_ukf(:,:,k+1)] = ukf_input(x_ukf(:,k), f_measured(k),in(:,k), P_ukf(:,:,k), f, h, Q, R); % Calls UKF function
end


tiledlayout(2,4);
nexttile;
plot(t, x_ukf(1,:)); hold all; plot(t,z);title("z"); legend("UKF","True"); xlabel("Time(s)"); ylabel("x");
nexttile;
plot(t, x_ukf(2,:)); hold all; plot(t,c0);title("c0"); legend("UKF","True");xlabel("Time(s)"); ylabel("v");
nexttile;
plot(t, x_ukf(3,:)); hold all; plot(t,k0);title("k0"); legend("UKF","True");
nexttile;
plot(t, x_ukf(4,:)); hold all; plot(t,alpha);title("alpha"); legend("UKF","True"); xlabel("Time(s)"); ylabel("x");
nexttile;
plot(t, x_ukf(5,:)); hold all; plot(t,beta);title("beta"); legend("UKF","True");xlabel("Time(s)"); ylabel("v");
nexttile;
plot(t, x_ukf(6,:)); hold all; plot(t,gamma);title("gamma"); legend("UKF","True");
nexttile;
plot(t, x_ukf(7,:)); hold all; plot(t,A);title("A"); legend("UKF","True"); xlabel("Time(s)"); ylabel("x");
nexttile;
plot(t, x_ukf(8,:)); hold all; plot(t,x0);title("x0"); legend("UKF","True"); xlabel("Time(s)"); ylabel("x");
