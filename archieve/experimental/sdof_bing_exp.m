clc;
close all;

T = 0.0147;

data = readtable(['data_voltagevarying_0volt.csv']);

t = data(:,2);
t = t{:,:};
f_measured = data(:,6);
f_measured= f_measured{:,:};

x0 = data(:,3);
x0 = x0{:,:};
x1 = data(:,4);
x1 = x1{:,:};

order = 3;
framelen = 11;
x1 = sgolayfilt(x1,order,framelen);

x_True = [x0 x1]';

%% DEFINING NON-LINEAR EQNS %%%
n = 2; 
m = 1;

%% PARAMETERS (FROM LLS) %%
cd = 352.5194;    
fc = 5.5396;
f0 = -8.2707;
mass = 0.1

f = @(x,u) [
    x(1)+x(2)*T;
    x(2)+(cd*x(2)+fc*sign(x(2))+f0)*T/mass
    ];

h = @(x,u)[
    cd*x(2)+fc*sign(x(2))+f0
];

Q = 1e-3*eye(n); 
R = 1e-3*eye(m); 
P_ukf = 10*Q;
f_ukf = zeros(m, length(t)); 
x_ukf = zeros(n, length(t)); 

%% UNSCENTED KALMAN FILTER%%
for k = 1:length(t)-1
    f_ukf(k+1) = h(x_ukf(:,k),0);
    [x_ukf(:,k+1), P_ukf(:,:,k+1)] = ukf_input(x_ukf(:,k), f_measured(k),0, P_ukf(:,:,k), f, h, Q, R); 
end

f_true = zeros(m, length(t));
for k = 1:length(t)-1
    f_true(k+1) = cd*x1(k)+fc*sign(x1(k))+f0;
end

tiledlayout(2,2);
nexttile;
plot(t, x_ukf(1,:)); hold all; plot(t,x_True(1,:));title("Position"); legend("UKF","True");
nexttile;
plot(t, x_ukf(2,:)); hold all; plot(t,x_True(2,:));title("Velocity"); legend("UKF","True");
% nexttile;
% plot(x_ukf(2,:),f_ukf); hold all; plot(x_True(2,:),f_measured);title("FV"); legend("UKF","True");
nexttile;
plot(t,f_true); hold all; ;plot(t,f_ukf); plot(t,f_measured);title("FV"); legend("Bingham Fitting (LLS)","Bingham Fitting (UKF)", "Measured");

