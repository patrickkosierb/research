%% ur5 general momentum state space kf %%
% init
% run_model

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
q_noisey = noisey(q_true,m,Rq,t_len); 
qd_noisey = noisey(qd_true,m,Rqd,t_len);
Q = 10e-6*eye(N); 
Q_tau = 20e-2; %3 works too
Q(7,7) = Q_tau;
Q(8,8) = Q_tau;
Q(9,9) = Q_tau;
Q(10,10) = Q_tau;
Q(11,11) = Q_tau;
Q(12,12) = Q_tau;
P = zeros(N,N,t_len);
P(:,:,1) = Q;
% state meas and input
p_meas = zeros(dof,t_len);
for k=1:t_len
    p_meas(:,k) = massMatrix(ur5_rbt,q_noisey(:,k))*qd_noisey(:,k);
end
x_estgm = zeros(dof*n,t_len);
u = zeros(dof*n,t_len);

%% ss
A = [zero_dxd I_dxd; zero_dxd zero_dxd];
B = [I_dxd zero_dxd;zero_dxd zero_dxd];
C = [I_dxd zero_dxd];

% discretize zoh
A = eye(dof*n)+A*Ts;
B = B*Ts;

%% kf
tic;
for k=1:t_len-1
    u(1:dof,k) = tau_jint_out(:,k)+ur5.coriolis(q_noisey(:,k)',qd_noisey(:,k)')'*qd_noisey(:,k)-gravityTorque(ur5_rbt,q_noisey(:,k));
    x_pred=A*x_estgm(:,k)+B*u(:,k);
    P_pred = A*P(:,:,k)*A'+Q;
    K = P_pred*C'/(C*P_pred*C' + Rqd);
    x_estgm(:,k+1) = x_pred+K*(p_meas(:,k)-C*x_pred);
    P(:,:,k+1) = (eye(dof*n) - K*C)*P_pred*(eye(dof*n) - K*C)' + K*Rqd*K'; 
end
GMKF_RT = toc;

errortau_gmkf = tau_ext_ideal(6,:)+x_estgm(12,:);
rmsetau_gmkf = sqrt(mean((errortau_gmkf).^2));

