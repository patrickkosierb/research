%% ur5_sim with fomo %%
% init
% run_model

%% first-order momentum observer %%

% import ur5 from RTB(Peter Corke) for Coriolis  
mdl_ur5
for i = 1:6
    ur5.links(i).Jm = 0;  % Example motor inertia for all links
    ur5.links(i).G = -9.81;
end
ur5.links(1).I = [0.010267495893, 0.010267495893, 0.00666, 0,0,0];
ur5.links(2).I = [0.22689067591, 0.22689067591, 0.0151074, 0,0,0];
ur5.links(3).I = [0.049443313556, 0.049443313556, 0.004095, 0,0,0];
ur5.links(4).I = [0.111172755531, 0.111172755531, 0.21942, 0,0,0];
ur5.links(5).I = [0.111172755531, 0.111172755531, 0.21942, 0,0,0];
ur5.links(6).I = [0.0171364731454, 0.0171364731454, 0.033822, 0,0,0];

% add noise to position & velcocity
m=6;
Rq = 10e-9*eye(m); 
Rqd = 10e-5*eye(m); 
q_noisey = noisey(q_true,m,Rq,t_len);
qd_noisey = noisey(qd_true,m,Rqd,t_len);

% gmo param
r = zeros(6,t_len);
p0 = massMatrix(ur5_rbt,q_noisey(:,1))*qd_noisey(:,1);
K=13; %2.6
tau_c = 1/K;
int_running = 0;

% fomo loop
tic;
for k=1:t_len-1
    integral = tau_jint_out(:,k+1)+ur5.coriolis(q_noisey(:,k+1)',qd_noisey(:,k+1)')'*qd_noisey(:,k+1)-gravityTorque(ur5_rbt,q_noisey(:,k+1))+r(:,k);
    int_running = int_running+integral.*Ts;
    r(:,k+1)=K*(massMatrix(ur5_rbt,q_noisey(:,k+1))*qd_noisey(:,k+1)-int_running-p0);
end
FOMO_RT = toc;

% figure;
% plot(t,tau_ext_ideal(6,:));
% hold on;
% plot(t,-r(6,:));
% legend("Actual", "FOMO")

errortau_fomo = tau_ext_ideal(6,:)+r(6,:);
rmsetau_fomo = sqrt(mean((errortau_fomo).^2));
