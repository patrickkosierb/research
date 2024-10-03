clc;close all;clear;

addpath("C:\Users\pk\Documents\code\research\data");
addpath("C:\Users\pk\Documents\code\research\general");

% pre-processing
data = readtable('data_voltagevarying.csv');
t = data(:,1);
t = t{:,:};
x = data(:,2);
x = x{:,:};
xdot = data(:,3);
xdot = xdot{:,:};
voltage = data(:,4);
voltage = voltage{:,:};
f_measured = data(:,5);
f_measured= f_measured{:,:};

% filtering noisey velocity
order = 3;
framelen = 11;
xdot = sgolayfilt(xdot,order,framelen);

% plotting data 
% figure;
% tiledlayout(4,1);
% nexttile;
% plot(t, x); hold all; title("position"); xlabel("time(s)"); ylabel("m");
% nexttile;
% plot(t, xdot); hold all; title("velocity"); xlabel("time(s)"); ylabel("m/s");
% nexttile;
% plot(t, voltage); hold all; title("voltage");xlabel("time(s)"); ylabel("V");
% nexttile;
% plot(t, f_measured); hold all; title("force");xlabel("time(s)"); ylabel("V");

% num. of states & measurement + sampling time
n = 2;
m = 1;
T = 0.0147;

% 'true' parameters from nlls
n_po = 2;
c = 5735.324123113438;
k = 222.1732327530771;
alpha = 134.46098993419324;
beta = 1519.1475673855064;
gamma = 1485.0310631782363;
A = -0.15334767682272726;
x0 = -0.0029967590221754714;

% state func.
f = @(x,u)[
    x(1)+(-gamma*abs(u(2))*x(1)*abs(x(1))^(n_po-1)-beta*u(2)*abs(x(1))^n_po+A*u(2))*T % z
    x(2); % c0
];
% measurement func.
h = @(x,u)[
    x(2)*u(2)+k*(u(1)-x0)+alpha*x(1); % force
];

n_iter = round(length(t)/3);
left_bound = 1;

Q = 1e-7*eye(n); 
R = 1e-4*eye(m); 
P_ukf = Q;

x_est = zeros(n, n_iter-left_bound-1);
x_est(:,1) = [0.1,1];
f_est = zeros(m, n_iter-left_bound-1);
in = [x,xdot]';

% ukf loop
for k = 1:n_iter-1
    [x_est(:,k+1), P_ukf] = ukf(x_est(:,k), f_measured(k), in(:,k), P_ukf, f, h, Q, R);
    f_est(k+1) = h(x_est(:,k+1),in(:,k));
end

% round(length(t)/5);

plot(t(left_bound:n_iter),x_est(1,left_bound:n_iter));
title("param c")
figure;
plot(t(left_bound:n_iter),f_measured(left_bound:n_iter));
hold on;
plot(t(left_bound:n_iter),f_est(left_bound:n_iter));
title("force")





