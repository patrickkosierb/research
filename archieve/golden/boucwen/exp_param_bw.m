clc;close all;clear

T = 0.0147;
tic;

data = readtable(['v250000_a50_121mv2_clean.csv']);

t = data(:,2);
t = t{:,:};

f_measured = data(:,6);
f_measured= f_measured{:,:};

x = data(:,3);
x = x{:,:}*100;
xdot = data(:,4);
xdot = xdot{:,:};

order = 3;
framelen = 11;
xdot = sgolayfilt(xdot,order,framelen)*100;

in = [x xdot]';
n = 3;
m = 1;

%% PARAMETERS %%
n_po = 2;

c0= 8.86846727390849;
k0= 0.597632140788303;
alpha= 30.837990991911244;
beta= -954.4497438208226;
gamma= 972.7733104878047;
A= 7.622000649825172;
x0= 5.388211694664788;

c0= ones(1,length(t))*c0;
k0= ones(1,length(t))*k0;
z = zeros(1,length(t));

x_True = [z' c0' k0']';

f = @(x,u)[
    x(1)-(-gamma*abs(u(2))*x(1)*abs(x(1))^(n_po-1)-beta*u(2)*abs(x(1))^n_po+A*u(2))*T
    x(2); %c0
    x(3); %k0
];

% MEASUREMENT VECTOR
h = @(x,u)[
    x(2)*u(2)+x(3)*(u(1)-x0)+alpha*x(1);
];

% NOISE
Q = 1e-8*eye(n); 
R = 1e-4*eye(m); 
P_ukf = 10*Q;

x_ukf = zeros(n,length(t));
x_ukf(:,1)=[0.1,1,1];

%% UNSCENTED KALMAN FILTER%%
for k = 1:length(t)-1
    [x_ukf(:,k+1), P_ukf(:,:,k+1)] = ukf_input(x_ukf(:,k), f_measured(k),in(:,k), P_ukf(:,:,k), f, h, Q, R); % Calls UKF function
end

f_lls = zeros(1,length(t));
f_ukf = zeros(1,length(t));
z_lls = zeros(1,length(t));
c0_lls= 8.86846727390849;
k0_lls= 0.597632140788303;

c0_ukf = mean(x_ukf(2,1000:length(t)))
k0_ukf = mean(x_ukf(3,1000:length(t)))

disp(toc)

z_lls(1) =  0.1;
for k = 1:length(t)-1
    z_lls(k+1) = z_lls(k)+(-gamma*abs(in(2,k))*z_lls(k)*abs(z_lls(k))^(n_po-1)-beta*in(2,k)*abs(z_lls(k))^n_po+A*in(2,k))*T;
    f_lls(k) = c0_lls*in(2,k)+k0_lls*(in(1,k)-x0)+alpha*z_lls(k);
    f_ukf(k) = c0_ukf*in(2,k)+k0_ukf*(in(1,k)-x0)+alpha*x_ukf(1,k);
end

squaredDifferences = (f_measured' - f_ukf).^2;
mse = mean(squaredDifferences);
rmse_ukf = sqrt(mse);

squaredDifferences = (f_measured' - f_lls).^2;
mse = mean(squaredDifferences);
rmse_lls = sqrt(mse);

paramName = {'c0', 'k0'};
paramUKF = [c0_ukf, k0_ukf];
paramLLS = [c0_lls, k0_lls];
Tparam = table(paramName', paramUKF', paramLLS', 'VariableNames', {'Parameter', 'UKF', 'LLS'});
disp(Tparam);

method = {'UKF', 'LLS'};
rmseCol = [rmse_ukf, rmse_lls];
Trmse = table(method', rmseCol', 'VariableNames', {'Method', 'RMSE'});
disp(Trmse);

tiledlayout(2,2);
nexttile;
plot(t, x_ukf(1,:)); hold all;plot(t, z_lls); title("Evolutionary Variable Z"); legend("UKF", "LLS"); xlabel("Time(s)"); ylabel("Ns/cm");
nexttile;
plot(t, x_ukf(2,:)); hold all; plot(t,ones(1,length(t))*c0_ukf);  plot(t,ones(1,length(t))*c0_lls);title("c0 Parameter ID"); legend("UKF","c0 mean", "LLS"); xlabel("Time(s)"); ylabel("Ns/cm");
nexttile;
plot(t, x_ukf(3,:)); hold all; plot(t,ones(1,length(t))*k0_ukf);  plot(t,ones(1,length(t))*k0_lls);title("k0 Parameter ID"); legend("UKF","k0 mean", "LLS"); xlabel("Time(s)"); ylabel("N/cm");
nexttile;
plot(t, f_ukf, 	'Color',[0 0.4470 0.7410]); hold all;  plot(t,f_lls,'Color',[0.9290 0.6940 0.1250]); plot(t,f_measured,'r'); title("Measured vs. UKF vs. LLS"); legend("UKF","LLS","True"); xlabel("Time(s)"); ylabel("N");






