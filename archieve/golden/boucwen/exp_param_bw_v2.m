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
n = 8;
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
alpha= ones(1,length(t))*alpha;
gamma= ones(1,length(t))*gamma;
beta= ones(1,length(t))*beta;
A= ones(1,length(t))*A;
x0= ones(1,length(t))*x0;


z = zeros(1,length(t));

x_True = [z' c0' k0' alpha' gamma' beta' A' x0']';

f = @(x,u)[
    x(1)-(-x(5)*abs(u(2))*x(1)*abs(x(1))^(n_po-1)-x(6)*u(2)*abs(x(1))^n_po+x(7)*u(2))*T
    x(2); %c0
    x(3); %k0
    x(4);
    x(5);
    x(6);
    x(7);
    x(8);
];

% MEASUREMENT VECTOR
h = @(x,u)[
    x(2)*u(2)+x(3)*(u(1)-x(8))+x(4)*x(1);
];

% NOISE
Q = 1e-8*eye(n); 
R = 1e-4*eye(m); 
P_ukf = 10*Q;

x_ukf = zeros(n,length(t));
x_ukf(:,1)=[0.1,1,1,10,0,0,0, 0];

%% UNSCENTED KALMAN FILTER%%
for k = 1:length(t)-1
    [x_ukf(:,k+1), P_ukf(:,:,k+1)] = ukf_input(x_ukf(:,k), f_measured(k),in(:,k), P_ukf(:,:,k), f, h, Q, R); % Calls UKF function
end

f_lls = zeros(1,length(t));
f_ukf = zeros(1,length(t));
z_lls = zeros(1,length(t));
c0_lls= 8.86846727390849;
k0_lls= 0.597632140788303;
alpha_lls = 30.837990991911244;
gamma_lls= 972.7733104878047;
beta_lls= -954.4497438208226;
A_lls= 7.622000649825172;
x0_lls= 5.388211694664788;


c0_ukf = mean(x_ukf(2,1000:length(t)))
k0_ukf = mean(x_ukf(3,1000:length(t)))
alpha_ukf = mean(x_ukf(4,1000:length(t)))
gamma_ukf = mean(x_ukf(5,1000:length(t)))
beta_ukf = mean(x_ukf(6,1000:length(t)))
A_ukf = mean(x_ukf(7,1000:length(t)))
x0_ukf = mean(x_ukf(8,1000:length(t)))


disp(toc)

z_lls(1) =  0.1;
for k = 1:length(t)-1
    z_lls(k+1) = z_lls(k)+(-gamma_lls*abs(in(2,k))*z_lls(k)*abs(z_lls(k))^(n_po-1)-beta_lls*in(2,k)*abs(z_lls(k))^n_po+A_lls*in(2,k))*T;
    f_lls(k) = c0_lls*in(2,k)+k0_lls*(in(1,k)-x0_lls)+alpha_lls*z_lls(k);
    f_ukf(k) = c0_ukf*in(2,k)+k0_ukf*(in(1,k)-x0_ukf)+alpha_ukf*x_ukf(1,k);
end

squaredDifferences = (f_measured(3000:length(t))' - f_ukf(3000:length(t))).^2;
mse = mean(squaredDifferences);
rmse_ukf = sqrt(mse);

squaredDifferences = (f_measured(3000:length(t))' - f_lls(3000:length(t))).^2;
mse = mean(squaredDifferences);
rmse_lls = sqrt(mse);

paramName = {'c0', 'k0', 'alpha', 'gamma', 'beta', 'A', 'x0'};
paramUKF = [c0_ukf, k0_ukf, alpha_ukf gamma_ukf beta_ukf A_ukf x0_ukf];
paramLLS = [c0_lls, k0_lls, alpha_lls gamma_lls beta_lls A_lls x0_lls];
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
plot(t, x_ukf(4,:)); hold all; plot(t,ones(1,length(t))*alpha_ukf);  plot(t,ones(1,length(t))*alpha_lls);title("alpha Parameter ID"); legend("UKF","alpha mean", "LLS"); xlabel("Time(s)"); ylabel("N/cm");
figure;
tiledlayout(2,2);
nexttile;
plot(t, x_ukf(5,:)); hold all; plot(t,ones(1,length(t))*gamma_ukf);  plot(t,ones(1,length(t))*gamma_lls);title("gamma Parameter ID"); legend("UKF","gamma mean", "LLS"); xlabel("Time(s)"); ylabel("N/cm");
nexttile;
plot(t, x_ukf(6,:)); hold all; plot(t,ones(1,length(t))*beta_ukf);  plot(t,ones(1,length(t))*beta_lls);title("beta Parameter ID"); legend("UKF","beta mean", "LLS"); xlabel("Time(s)"); ylabel("N/cm");
nexttile;
plot(t, x_ukf(7,:)); hold all; plot(t,ones(1,length(t))*A_ukf);  plot(t,ones(1,length(t))*A_lls);title("A Parameter ID"); legend("UKF","A mean", "LLS"); xlabel("Time(s)"); ylabel("N/cm");
nexttile;
plot(t, x_ukf(8,:)); hold all; plot(t,ones(1,length(t))*x0_ukf);  plot(t,ones(1,length(t))*x0_lls);title("x0 Parameter ID"); legend("UKF","x0 mean", "LLS"); xlabel("Time(s)"); ylabel("N/cm");


figure;
plot(t, f_ukf, 	'Color',[0 0.4470 0.7410]); hold all;  plot(t,f_lls,'Color',[0.9290 0.6940 0.1250]); plot(t,f_measured,'r'); title("Measured vs. UKF vs. LLS"); legend("UKF","LLS","True"); xlabel("Time(s)"); ylabel("N");

