clc;close all;clear

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
figure;
tiledlayout(4,1);
nexttile;
plot(t, x); hold all; title("position"); xlabel("time(s)"); ylabel("m");
nexttile;
plot(t, xdot); hold all; title("velocity"); xlabel("time(s)"); ylabel("m/s");
nexttile;
plot(t, voltage); hold all; title("voltage");xlabel("time(s)"); ylabel("V");
nexttile;
plot(t, f_measured); hold all; title("force");xlabel("time(s)"); ylabel("V");

% num. of states & measurement + sampling time
n = 4;
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

% lls force estimate
f_lls = zeros(size(t));
z_lls = zeros(size(t));
z_lls(1) =  0.1;
for k = 1:length(t)-1
    z_lls(k+1) = z_lls(k)+(-gamma*abs(xdot(k))*z_lls(k)*abs(z_lls(k))^(n_po-1)-beta*xdot(k)*abs(z_lls(k))^n_po+A*xdot(k))*T;
    f_lls(k) = c*xdot(k)+k*(xdot(k)-x0)+alpha*z_lls(k);
end
figure;
plot(t,f_lls);

% state func.
f = @(x,u)[
    x(1)+(-gamma*abs(u(2))*x(1)*abs(x(1))^(n_po-1)-beta*u(2)*abs(x(1))^n_po+A*u(2))*T % z
    x(2); % c0
    x(3); % k0
    x(4); % alpha
];
% measurement func.
h = @(x,u)[
    x(2)*u(2)+x(3)*(u(1)-x0)+x(4)*x(1); % force
];

Q = 1e-7*eye(n); 
R = 1e-4*eye(m); 
P_ukf = 10*Q;

x_est = zeros(n, length(t));
x_est(:,1) = [0.1,1,1,1];

in = [x;xdot];

% ukf
for k = 1:length(t)-1
    [x_est(:,k+1), P_ukf(:,:,k+1)] = ukf_input(x_est(:,k), f_measured(k), in(:,k), P_ukf(:,:,k), f, h, Q, R); 
end