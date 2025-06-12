%% init virtual sensor sim %%
% clear; clc; close all;

addpath(pwd+"/helpers/")
addpath(pwd+"/Geometry/")

% simulation param 
ur5e = importrobot("ur5_block",DataFormat="column")
runtime = 2;
Ts = 0.01;
t = 0:Ts:runtime;
t_len = length(t);
TIM_MAG = 1/Ts;
% ext torque 
delay = 1*TIM_MAG;
imp = 0.5;
tau_ext = 0.5;  %0.5
% traj. param 
waypoints = [0 0 0 0 0 0;0 0 0 0 0 0]';
time = [0,1];
v_bound = zeros(6,runtime);
fext = zeros(6,10); 

mdl_ur5

for i = 1:6
    ur5.links(i).Jm = 0;  % Example motor inertia for all links
    ur5.links(i).G = 9.81;
end

ur5.links(1).I = [0.010267495893, 0.010267495893, 0.00666, 0,0,0];
ur5.links(2).I = [0.22689067591, 0.22689067591, 0.0151074, 0,0,0];
ur5.links(3).I = [0.049443313556, 0.049443313556, 0.004095, 0,0,0];
ur5.links(4).I = [0.111172755531, 0.111172755531, 0.21942, 0,0,0];
ur5.links(5).I = [0.111172755531, 0.111172755531, 0.21942, 0,0,0];
ur5.links(6).I = [0.0171364731454, 0.0171364731454, 0.033822 0,0,0];

out = sim("ur5_sim.slx",runtime);

% process output from sim
q_true = out.q_out.Data';
qd_true = out.qdot_out.Data';
qddot_true = out.qddot_out.Data';
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
    C = ur5.coriolis(q_true(:,k)',qd_true(:,k)');
    Cv = C*qd_true(:,k);
    tau_ext_ideal(:,k) = tau_jint_out(:,k)-tau_G_out(k,:)'-tau_C_out(k,:)'-tau_M_out(:,:,k)*qddot_true(:,k);
    tau_ext_ideal_rtb(:,k) = tau_jint_out(:,k)-tau_G_out(k,:)'-Cv-tau_M_out(:,:,k)*qddot_true(:,k);
end

plot(t,tau_ext_ideal(6,:))
hold on;
plot(t,tau_ext_ideal_rtb(6,:))
legend("rstb","rtb")

error = tau_ext_ideal(6,:)-tau_ext_ideal_rtb(6,:);
rmse = sqrt(mean((error).^2))