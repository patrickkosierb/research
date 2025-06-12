%% ur5 control with torque from urSIM data %%
init
%% run model
out = sim(sim_to_run,RUNTIME);

% process output from sim
q_true = out.q_out.Data';
qd_true = out.qdot_out.Data';
qdd_true = out.qddot_out.Data';
tau_jint_out = out.tau_jint_out.Data';

% numemrical approx of accel
% qdd_calc = zeros(6,t_span);
% for i=1:t_span-1
%     qdd_calc(:,i) = (qd_actual(:,i+1)-qd_actual(:,i))/Ts;
% end
% % smooth accel
% qdd_calc = sgolayfilt(qdd_calc, 3, 5);

% calculate external torque
tau_ext_ideal = zeros(6,t_span);
for k=1:t_span
    tau_ext_ideal(:,k) = tau_jint_out(:,k)-gravityTorque(ur5_rbt,q_true(:,k))-velocityProduct(ur5_rbt,q_true(:,k),qd_true(:,k))-massMatrix(ur5_rbt,q_true(:,k))*qdd_true(:,k);
end

% plot(time,q_true)
% legend("1","2","3","4","5","6");
plot(time,tau_ext_ideal)
legend("1","2","3","4","5","6");
% plot(time,tau_ext_ideal)
% figure;
% plot(time,tau_jint_out)