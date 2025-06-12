%% ur5 general momentum state space kf %%
run_model

%% general param
dof = 6;
n = 2;
m=6;
N = dof*n;
zero_dxd = zeros(dof,dof);
I_dxd = eye(dof,dof);

%% kf param
% add noise sim measurements
Rq = 10e-9*eye(m); 
Rqd = 10e-5*eye(m);
q_noisey = noisey(q_true,m,Rq,t_span); 
qd_noisey = noisey(qd_true,m,Rqd,t_span);
Q = 10e-6*eye(N); 
Q_tau = 20e-2; %3 works too
Q(7,7) = Q_tau;
Q(8,8) = Q_tau;
Q(9,9) = Q_tau;
Q(10,10) = Q_tau;
Q(11,11) = Q_tau;
Q(12,12) = Q_tau;
P = zeros(N,N,t_span);
P(:,:,1) = Q;
% state meas and input
p_meas = zeros(dof,t_span);
for k=1:t_span
    p_meas(:,k) = massMatrix(ur5_rbt,q_noisey(:,k))*qd_noisey(:,k);
end
x_est = zeros(dof*n,t_span);
u = zeros(dof*n,t_span);

%% ss
A = [zero_dxd I_dxd; zero_dxd zero_dxd];
B = [I_dxd zero_dxd;zero_dxd zero_dxd];
C = [I_dxd zero_dxd];
% discretize zoh
A = eye(dof*n)+A*Ts;
B = B*Ts;

%% kf
tic;
for k=1:t_span-1
    u(1:dof,k) = tau_jint_out(:,k)+ur5e.coriolis(q_noisey(:,k)',qd_noisey(:,k)')'*qd_noisey(:,k)-gravityTorque(ur5_rbt,q_noisey(:,k));
    x_pred=A*x_est(:,k)+B*u(:,k);
    P_pred = A*P(:,:,k)*A'+Q;
    K = P_pred*C'/(C*P_pred*C' + Rqd);
    x_est(:,k+1) = x_pred+K*(p_meas(:,k)-C*x_pred);
    P(:,:,k+1) = (eye(dof*n) - K*C)*P_pred*(eye(dof*n) - K*C)' + K*Rqd*K'; 
end
toc

figure;
plot(time,tau_ext_ideal(6,:));
hold on;
plot(time,-x_est(12,:))
legend("Actual", "GMKF")

errortau_gmkf = tau_ext_ideal(6,:)+x_est(12,:);
rmsetau_gmkf = sqrt(mean((errortau_gmkf).^2))