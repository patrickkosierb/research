clc;clear;close all;

T=0.1;
t = 0:T:500;

n = 2; 
m = 1;

%% TRUE DYNAMICS %%
period = 15;
w = 2*pi*1/period;
x = (-0.0375*sin(w*t-pi/2)-0.0375)*100;
xdot = -0.0375*w*cos(w*t-pi/2)*100;
xddot = 0.0375*w^2*sin(w*t-pi/2)*100;



%% PARAMETERS %%
% c0= 11.8;
% k0= 0.6;
% alpha= 44.5;
% beta= -3062;
% gamma= 3078;
% A= 2.8;
% x0= 7.2;
% n_po = 2;

c0= 8.86846727390849;
k0= 0.597632140788303;
alpha= 30.837990991911244;
beta= -954.4497438208226;
gamma= 972.7733104878047;
A= 7.622000649825172;
x0= 5.388211694664788;
n_po = 2;

z = zeros(1,length(t)-1);
z(1)=1;
for k = 1:length(t)-1
    z(k+1) = z(k)+(-gamma*abs(xdot(k))*z(k)*abs(z(k))^(n_po-1)-beta*xdot(k)*abs(z(k))^n_po+A*xdot(k))*T

end
z;
in = [x' xdot' z']';

c0= ones(1,length(t))*c0;
k0= ones(1,length(t))*k0;

x_True = [c0' k0']';

f = @(x,u)[
    % x(1)+(-gamma*abs(u(2))*x(1)*abs(x(1))^(n_po-1)-beta*u(2)*abs(x(1))^n_po+A*u(2))*T
    x(1); %c0
    x(2); %k0
];

% MEASUREMENT VECTOR
h = @(x,u)[
    x(1)*u(2)+x(2)*(u(1)-x0)+alpha*u(3);
];

% NOISE
Q = 1e-2*eye(n); 
R = 1e-2*eye(m); 
P_ukf = 10*Q;

rng('shuffle')
random_seed = rng;
w = mvnrnd(zeros(length(t),n),Q)'; 
v = mvnrnd(zeros(length(t),m),R)'; 

f_measured = zeros(m,length(t));
x_ukf = zeros(n,length(t));


%% UNSCENTED KALMAN FILTER%%
for k = 1:length(t)-1
    x_True(1,k) = x_ukf(1,k);
    f_measured(k+1) = h(x_True(:,k),in(:,k))+v(k);
    [x_ukf(:,k+1), P_ukf(:,:,k+1)] = ukf_input(x_ukf(:,k), f_measured(k),in(:,k), P_ukf(:,:,k), f, h, Q, R); % Calls UKF function
end

c0_ukf = mean(x_ukf(2,500:length(t)))
k0_ukf = mean(x_ukf(3,500:length(t)))
f_ukf = zeros(1,length(t));

for k = 1:length(t)-1
   f_ukf(k)= c0_ukf*in(2,k)+k0_ukf*(in(1,k)-x0)+alpha*x_ukf(1,k);
end


tiledlayout(2,2);
nexttile;
plot(t, x_ukf(1,:)); hold all; title("Evolutionary Variable"); ;xlabel("Time(s)"); ylabel("");
nexttile;
plot(t, x_ukf(2,:)); hold all; plot(t, ones(1,length(t))*c0_ukf); plot(t,x_True(2,:));title("c0 Parameter ID"); legend("UKF","c0 mean","True"); xlabel("Time(s)"); ylabel("Ns/cm");
nexttile;
plot(t, x_ukf(3,:)); hold all; plot(t, ones(1,length(t))*k0_ukf); plot(t,x_True(3,:));title("k0 Parameter ID"); legend("UKF","k0 mean","True");xlabel("Time(s)"); ylabel("N/cm");
nexttile;
plot(t, f_ukf); hold all; plot(t, f_measured); title("Simulated Measured vs. UKF");legend("UKF","True"); xlabel("Time(s)"); ylabel("N");







