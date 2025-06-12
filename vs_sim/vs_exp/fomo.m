%% ur5_dyanmics with fomo %%
run_model

%% first-order momentum observer %%

% add noise to position & velcocity
m=6;
Rq = 10e-9*eye(m); 
Rqd = 10e-5*eye(m); 
q_noisey = noisey(q_true,m,Rq,t_span);
qd_noisey = noisey(qd_true,m,Rqd,t_span);

% gmo param
r = zeros(6,t_span);
p0 = massMatrix(ur5_rbt,q_noisey(:,1))*qd_noisey(:,1);
K=13; %2.6
tau_c = 1/K;
int_running = 0;

% fomo loop
tic;
for k=1:t_span-1
    integral = tau_jint_out(:,k+1)+ur5e.coriolis(q_noisey(:,k+1)',qd_noisey(:,k+1)')'*qd_noisey(:,k+1)-gravityTorque(ur5_rbt,q_noisey(:,k+1))+r(:,k);
    int_running = int_running+integral.*Ts;
    r(:,k+1)=K*(massMatrix(ur5_rbt,q_noisey(:,k+1))*qd_noisey(:,k+1)-int_running-p0);
end
toc

figure;
plot(time,tau_ext_ideal(6,:));
hold on;
plot(time,-r(6,:));
legend("Actual", "FOMO")

errortau_fomo = tau_ext_ideal(6,:)+r(6,:);
rmsetau_fomo = sqrt(mean((errortau_fomo).^2))
