%% init sim %%
clear; clc; close all;
addpath(pwd+"/helpers/")
addpath(pwd+"/Geometry/")
addpath(pwd+"/vs_exp/")

% preprocessing urSIM data
sim_to_run = "ur5_dynamic_tauctrl_pid.slx";
ur5_rbt = importrobot("ur5_block",DataFormat="column");
mdl_ur5e
data = readtable('noload_100hz_enc.csv');
q_actual = [data.actual_q_0(:,1)  data.actual_q_1(:,1) data.actual_q_2(:,1) data.actual_q_3(:,1) data.actual_q_4(:,1) data.actual_q_5(:,1) ]';
qd_actual = [data.actual_qd_0(:,1)  data.actual_qd_1(:,1) data.actual_qd_2(:,1) data.actual_qd_3(:,1) data.actual_qd_4(:,1) data.actual_qd_5(:,1) ]';
time  = data.timestamp(:,1)';
t_span = length(time);
Ts=0.01;

% traj. for tau ctrl
bspline = [1,4];
% joint modes
mode_unlock = zeros(1,6);
mode_lock = ones(1,6);
RUNTIME = 4;
time= 0:Ts:4;
t_span = length(time);

% ext torque
TIM_MAG = 1/Ts;
delay = 0.5*TIM_MAG;
imp = 0.5;
tau_ext = 5;
fext = zeros(6,10); 
fext(:,10) = [0,0,0,0,0,0];

