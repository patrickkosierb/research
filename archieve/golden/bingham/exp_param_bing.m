clc;close all;clear

T = 0.0147;
tic;

data = readtable(['v250000_a50_121mv_clean.csv']);

t = data(:,2);
t = t{:,:};
f_measured = data(:,6);
f_measured= f_measured{:,:};

x = data(:,3);
x = x{:,:};
xdot = data(:,4);
xdot = xdot{:,:};

order = 3;
framelen = 11;
xdot = sgolayfilt(xdot,order,framelen);

n = 3; 
m = 1;

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

x_ukf = zeros(n,length(t));

%% UNSCENTED KALMAN FILTER%%
for k = 1:length(t)-1
    [x_ukf(:,k+1), P_ukf(:,:,k+1)] = ukf_input(x_ukf(:,k), f_measured(k),xdot(k), P_ukf(:,:,k), f, h, Q, R); % Calls UKF function
end

c0 = 16.2435;
fc = 14.5179;   
f0 = -6.3714;

coefficients = polyfit(t, x_ukf(1,:), 1);
% Extract slope and intercept from the coefficients
slope = coefficients(1);
c0_ukf = coefficients(2);
% c0_ukf = mean(x_ukf(1,:))
fc_ukf = mean(x_ukf(2,:));
f0_ukf = mean(x_ukf(3,:));

f_ukf = zeros(1,length(t));
f_lls = zeros(1, length(t));

for k = 1:length(t)-1
    f_ukf(k) = c0_ukf*100*xdot(k)+fc_ukf*sign(xdot(k))+f0_ukf;
    f_lls(k) = c0*100*xdot(k)+fc*sign(xdot(k))+f0;
end
toc
squaredDifferences = (f_measured' - f_ukf).^2;
mse = mean(squaredDifferences);
rmse_ukf = sqrt(mse);

squaredDifferences = (f_measured' - f_lls).^2;
mse = mean(squaredDifferences);
rmse_lls = sqrt(mse);

paramName = {'Fc', 'F0', 'Cd'};
paramUKF = [fc_ukf, f0_ukf, c0_ukf];
paramLLS = [fc, f0, c0];
Tparam = table(paramName', paramUKF', paramLLS', 'VariableNames', {'Parameter', 'UKF', 'LLS'});
disp(Tparam);

method = {'UKF', 'LLS'};
rmseCol = [rmse_ukf, rmse_lls];
Trmse = table(method', rmseCol', 'VariableNames', {'Method', 'RMSE'});
disp(Trmse);

tiledlayout(2,2)
nexttile
plot(t, x_ukf(1,:));
hold on 
plot(t,ones(1,length(t))*c0_ukf)
plot(t,ones(1,length(t))*c0)
legend("c0 UKF", "c0 LOBF", "c0 LLS")
title("c0 Parameter ID")
xlabel("Time(s)"); ylabel("Ns/m");
nexttile
plot(t, x_ukf(2,:));
hold on 
plot(t, ones(1,length(t))*fc_ukf);
plot(t,ones(1,length(t))*fc)
legend("fc UKF", "fc mean","fc LLS")
title("fc Parameter ID")
xlabel("Time(s)"); ylabel("N");
nexttile
plot(t, x_ukf(3,:));
hold on 
plot(t, ones(1,length(t))*f0_ukf);
plot(t,ones(1,length(t))*f0)
legend("f0 UKF", "f0 mean", "f0 LLS")
title("f0 Parameter ID")
xlabel("Time(s)"); ylabel("N");
nexttile
plot(t, f_ukf);
hold on;
plot(t, f_lls);
plot(t, f_measured);
legend("UKF", "LLS", "True")
title("Measured vs. UKF vs. LLS")
xlabel("Time(s)"); ylabel("N");

period = 15;
omega = 2*pi*1/period;
x_sim = -0.0375*sin(omega*t-pi/2)-0.0375;

figure;
plot(t, x);
hold all;
plot(t,x_sim);
legend("Experimental input position", "Simulation input position");
title("Input Position")
xlabel("Time(s)"); ylabel("m");