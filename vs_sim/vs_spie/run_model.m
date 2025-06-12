%% ur5_sim %%

%% run model
out = sim("ur5_sim.slx",runtime);

% process output from sim
q_true = out.q_out.Data';
qd_true = out.qdot_out.Data';
qdd_true = out.qddot_out.Data';
tau_jint_out = out.tau_jint_out.Data';
tau_jext_out = out.tau_jext_out.Data;
tau_G_out = out.tau_G.Data;
tau_C_out = out.tau_C.Data;
tau_M_out = out.tau_M.Data;

% total external torques relative to world frame
j1_ext = squeeze(tau_jext_out(1:3,:,:));
j2_ext = squeeze(tau_jext_out(4:6,:,:));
j3_ext = squeeze(tau_jext_out(7:9,:,:));
j4_ext = squeeze(tau_jext_out(10:12,:,:));
j5_ext = squeeze(tau_jext_out(13:15,:,:));
j6_ext = squeeze(tau_jext_out(16:18,:,:));

% calculate external torque
tau_ext_ideal = zeros(6,t_len);
for k=1:t_len
    tau_ext_ideal(:,k) = tau_jint_out(:,k)-tau_G_out(k,:)'-tau_C_out(k,:)'-tau_M_out(:,:,k)*qdd_true(:,k);
end