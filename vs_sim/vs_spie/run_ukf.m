%% ur5 with ukf %%
% init
% run_model

%% unscented kalman filter %%

% state transition function
f = @(x,u)[
 x(1)+x(7)*Ts; % q 
 x(2)+x(8)*Ts;
 x(3)+x(9)*Ts;
 x(4)+x(10)*Ts;
 x(5)+x(11)*Ts;
 x(6)+x(12)*Ts;
 x(7)+x(13)*Ts; % qd 
 x(8)+x(14)*Ts;            
 x(9)+x(15)*Ts;
 x(10)+x(16)*Ts;
 x(11)+x(17)*Ts;
 x(12)+x(18)*Ts;
 x(13); % qdd 
 x(14);
 x(15);
 x(16);
 x(17);
 x(18);
 x(19); % ext f/t 
 x(20);
 x(21);
 x(22);
 x(23);
 x(24);
];

% measurement function
h = @(x)[
 x(1); % q 
 x(2);
 x(3);
 x(4);
 x(5);
 x(6);
 x(7); % qd
 x(8);
 x(9);
 x(10);
 x(11);
 x(12);
];

% ukf params
jnts = 6;
n=4;
m=6;

Q = 10e-6*eye(n*jnts); %e-6
Q(19,19) = 20e-5;
Q(20,20) = 20e-5;
Q(21,21) = 20e-5;
Q(22,22) = 20e-5;
Q(23,23) = 20e-5;
Q(24,24) = 20e-5;
Rq = 10e-9*eye(m); 
Rqd = 10e-5*eye(m); %e-5
Z = zeros(m,m);
R = [Rq Z;Z Rqd];
kappa = 0.5;
P = Q;

% add noise to position
q_noisey = noisey(q_true,m,Rq,t_len);
qd_noisey = noisey(qd_true,m,Rqd,t_len);
meas = [q_noisey;qd_noisey];
% init est
x_est = zeros(n*jnts, t_len);

% ukf loop
tic;
for k = 1:t_len-1
    [x_est(:,k+1), P] = ukf_qdd(ur5_rbt,x_est(:,k), meas(:,k),tau_jint_out(:,k), P, f, h, Q, R,kappa);
end
UKF_RT = toc;

% true ee tau vs est tau
% figure;
% tiledlayout(3,1);
% nexttile;
% plot(t,tau_ext_ideal(6,:))
% hold on;
% plot(t,x_est(24,:))
% xlabel("time(s)")
% ylabel("torque(Nm)")
% legend("\tau_{ext} true", "\tau_{ext} est")
% title("tool z-axis external torque")
% 
% nexttile;
% plot(t,q_noisey(6,:))
% hold on;
% plot(t,x_est(6,:))
% xlabel("time(s)")
% ylabel("Position(rad)")
% legend("q true", "q est")
% title("tool z-axis position")
% 
% nexttile;
% plot(t,qd_noisey(6,:))
% hold on;
% plot(t,x_est(12,:))
% xlabel("time(s)")
% ylabel("Velocity(rad/s)")
% legend("qd true", "qd est")
% title("tool z-axis vel")
% 
% figure;
% plot(t,tau_ext_ideal(6,:))
% hold on;
% plot(t,x_est(24,:))
% xlabel("time(s)")
% ylabel("torque(Nm)")
% legend("\tau_{ext} true", "\tau_{ext} est")
% title("tool z-axis external torque")


%% evaluation metrics
% norm = abs(norm(tau_ext_ideal(6,:))-norm(x_est(24,:)))

errortau = tau_ext_ideal(6,:)-x_est(24,:);
rmsetau_ukf = sqrt(mean((errortau).^2));

% errorq = q_true(6,:)-x_est(6,:);
% rmseq = sqrt(mean((errorq).^2))
% 
% errorqd = qd_true(6,:)-x_est(12,:);
% rmseqd = sqrt(mean((errorqd).^2))

