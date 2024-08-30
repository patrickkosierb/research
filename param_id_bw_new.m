clc;close all;clear

% pre-processing
data = readtable(['data_voltagevarying.csv']);
t = data(:,1);
t = t{:,:};
f_measured = data(:,5);
f_measured= f_measured{:,:};
x = data(:,2);
x = x{:,:}*100;
xdot = data(:,3);
xdot = xdot{:,:};
voltage = data(:,4);
voltage= voltage{:,:};

% filtering noisey velocity
order = 3;
framelen = 11;
xdot = sgolayfilt(xdot,order,framelen)*100;

% num. of states & measurement + sampling time
n = 4;
m = 1;
T = 0.0147;

% true parameters from nlls
n_po = 2;
c0= 8.86846727390849;
k0= 0.597632140788303;
alpha= 30.837990991911244;
beta= -954.4497438208226;
gamma= 972.7733104878047;
A= 7.622000649825172;
x0= 5.388211694664788;

% to-be estimated 
z = zeros(1,length(t));
c0_list= ones(1,length(t))*c0;
k0_list = ones(1,length(t))*k0;
alpha_list = ones(1,length(t))*alpha;


% non-linear state function 
f = @(x,u)[
    x(1)-(-gamma*abs(u(2))*x(1)*abs(x(1))^(n_po-1)-beta*u(2)*abs(x(1))^n_po+A*u(2))*T % z
    x(2); % c0
    x(3); % k0
    x(4); % alpha
];

% measurement function 
h = @(x,u)[
    x(2)*u(2)+x(3)*(u(1)-x0)+x(4)*x(1);
];

% ukf inputs 
Q = 1e-7*eye(n); 
R = 1e-4*eye(m); 
P_ukf = 10*Q;

x_True = [z' c0_list' k0_list' alpha_list']';
x_ukf = zeros(n,length(t));
f_ukf = zeros(1,length(t));
in = [x xdot]';

%initial estimates 
x_ukf(:,1)=[0.1,1,1,10]; 

tic;
for k = 1:length(t)-1 % ukf loop
    [x_ukf(:,k+1), P_ukf(:,:,k+1)] = ukf_input(x_ukf(:,k), f_measured(k), in(:,k), P_ukf(:,:,k), f, h, Q, R); 
end

c0_ukf = mean(x_ukf(2,length(t)))
k0_ukf = mean(x_ukf(3,length(t)))
alpha_ukf = mean(x_ukf(4,length(t)))

toc;

for k = 1:length(t)-1
    f_ukf(k) = c0_ukf*in(2,k)+k0_ukf*(in(1,k)-x0)+alpha_ukf*x_ukf(1,k);
end

figure;
tiledlayout(3,1);
nexttile;
plot(t, in(1,:)); hold all; title("position"); xlabel("Time(s)"); ylabel("cm");
nexttile;
plot(t, in(2,:)); hold all; title("velocity"); xlabel("Time(s)"); ylabel("cm/s");
nexttile;
plot(t, voltage); hold all; title("applies voltage");xlabel("Time(s)"); ylabel("V");

figure;
tiledlayout(3,1);
nexttile;
plot(t, x_ukf(2,:)); hold all; plot(t,ones(1,length(t))*c0_ukf); title("c0 Parameter ID"); legend("UKF","c0 mean"); xlabel("Time(s)"); ylabel("Ns/cm");
nexttile;
plot(t, x_ukf(3,:)); hold all; plot(t,ones(1,length(t))*k0_ukf); title("k0 Parameter ID"); legend("UKF","k0 mean"); xlabel("Time(s)"); ylabel("N/cm");
nexttile;
plot(t, x_ukf(4,:)); hold all; plot(t,ones(1,length(t))*alpha_ukf); title("alpha Parameter ID"); legend("UKF","alpha mean"); xlabel("Time(s)"); ylabel("N/cm");

figure;
plot(t, f_ukf, 	'Color',[0 0.4470 0.7410]); hold all;  plot(t,f_measured,'r'); title("Measured vs. UKF"); legend("UKF","True"); xlabel("Time(s)"); ylabel("N");



