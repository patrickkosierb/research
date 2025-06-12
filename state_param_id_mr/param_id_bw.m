% Title: MR damper state + parameter estimation using the UKF
% Author: Patrick Kosierb // kosierbp@mcmaster.ca


% num. of states & measurement + sampling time
n = 6;
m = 3;
T = 0.0147;
t_len = length(t);

mass=0.43; 
Fg = mass*9.81;

% parameters from nlls
n_po = 2;

% medium speed
c_init= 8.868481050830804;
k_init=  0.5976260444329669;
alpha_init= 30.83878044539226;
beta= -954.5609515601979;
gamma= 972.8848513296159;
A= 7.621731293239106;
x0=  5.388361027175781;

% high speed
% c_init= 8.23755253142408;
% k_init= 0.8130750954878393;
% alpha_init= 40.98154769711848;
% beta= -1636.5555418048746;
% gamma= 1659.541803637622;
% A= 7.92904317184689;
% x0= 5.051566084440801;

% low speed
% c_init= 106.59349480511545;
% % c_init = 135.8883; % derived from UKF (mean from 1 run)
% k_init= 1.252438020628908;
% alpha_init= 8.584955025214565;
% beta= 0.04881631881513377;
% gamma= 0.2683090143904821;
% A= 0.510826704302544;
% x0= -1.2080818103750128;

% state func.
f = @(x,u)[
    x(1)+x(2)*T; %position
    x(2)+x(3)*T; %velcoity
    ((x(5)*x(2)+x(6)*(x(1)-x0)+alpha_init*x(4))-Fg)/mass %accel
    x(4)+(-gamma*abs(x(2))*x(4)*abs(x(4))^(n_po-1)-beta*x(2)*abs(x(4))^n_po+A*x(2))*T; % z
    x(5); %c
    x(6); %k
];
 
% measurement func.
h = @(x,u)[
    x(1); %position
    x(2); %velcoity
    x(5)*x(2)+x(6)*(x(1)-x0)+alpha_init*x(4); % force
]; 

% Q=  [1e-2 0 0 0 0 0; 0 1e-2 0 0 0 0; 0 0 1e-2 0 0 0; 0 0 0 1e-2 0 0; 0 0 0 0 8e-6 0; 0 0 0 0 0 8e-4 ]; %trial 1
% Q= [1e-2 0 0 0 0 0; 0 1e-1 0 0 0 0; 0 0 1e-2 0 0 0; 0 0 0 1e-2 0 0; 0 0 0 0 1e-12 0; 0 0 0 0 0 1e-12]; %trial2
Q = [1e-2 0 0 0 0 0; 0 5e-1 0 0 0 0; 0 0 1e-2 0 0 0; 0 0 0 1e-4 0 0; 0 0 0 0 1e-2 0; 0 0 0 0 0 1e-3]; %trial 3 


%consult the author for calculation of measurement covariance
R = [2.6694e-9 0 0;
    0 2.4708e-2 0; 
    0 0 0.0336]; 

%% UKF start
P_ukf = Q;
x_est = zeros(n, t_len);
x_est(:,1) = [0;0;0;0; c_init;k_init];
z_est = zeros(m, t_len);
in = 0;
z_meas = [x_meas,xdot,f_measured]';

tic;
% ukf loop
for k = 1:t_len-1
    [x_est(:,k+1), P_ukf] = ukf(x_est(:,k), z_meas(:,k), in, P_ukf, f, h, Q, R);
    z_est(:,k+1) = h(x_est(:,k+1),in);
end
PID_RT = toc;
rmse_x_pid = sqrt(sum((x_meas'-z_est(1,:)).^2)/t_len)
rmse_xdot_pid = sqrt(sum((xdot'-z_est(2,:)).^2)/t_len)

figure;
tiledlayout(2,1);
nexttile;
plot(t,x_meas);
hold on;
plot(t,x_est(1,:));
title("position")
xlabel("Time(ms)")
ylabel("Position(cm)")
legend("Measured", "Estimated")
nexttile;
plot(t,xdot);
hold on;
plot(t,x_est(2,:));
title("velocity")
xlabel("Time(ms)")
ylabel("Velocity(cm/s")
legend("Measured", "Estimated")

figure;
plot(f_measured,xdot);
hold on;
plot(z_est(3,:),x_est(2,:));
title("force velocity")
xlabel("Measured Force (N)")
ylabel("Velocity(cm/s")
legend("Measured", "Estimated")
% 
% figure;
% plot(t,x_est(5,:));
% c_ukf = sum(x_est(5,:))/t_len
% title("c estimate")
% xlabel("Time (Ms)")
% ylabel("C(Ns/cm)")
% 
% figure;
% plot(t,x_est(6,:));
% c_ukf = sum(x_est(6,:))/t_len
% title("k estimate")
% xlabel("Time (Ms)")
% ylabel("C(Ns/cm)")


% figure;
% plot(t,x_est(6,:));
% k_ukf = sum(x_est(6,:))/t_len
% title("k estimate")
% 
% figure;
% plot(t,x_est(7,:));
% alpha_ukf = sum(x_est(7,:))/t_len
% title("alpha estimate")
% figure;
% plot(t,f_measured);
% hold on;
% plot(t,z_est(3,:));
% title("force")
