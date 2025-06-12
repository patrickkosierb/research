%% init sim %%
clear; clc; close all;
addpath(pwd+"/helpers/")
addpath(pwd+"/Geometry/")
addpath(pwd+"/vs_exp/")

% preprocessing urSIM data
sim_to_run = "ur5_dynamic_tauctrl.slx";
ur5_rbt = importrobot("ur5_block",DataFormat="column");
mdl_ur5e
data = readtable('noload_100hz_enc.csv');
q_actual = [data.actual_q_0(:,1)  data.actual_q_1(:,1) data.actual_q_2(:,1) data.actual_q_3(:,1) data.actual_q_4(:,1) data.actual_q_5(:,1) ]';
qd_actual = [data.actual_qd_0(:,1)  data.actual_qd_1(:,1) data.actual_qd_2(:,1) data.actual_qd_3(:,1) data.actual_qd_4(:,1) data.actual_qd_5(:,1) ]';
time  = data.timestamp(:,1)';
t_span = length(time);
Ts=0.01;

if sim_to_run == "ur5_sim.slx"
    RUNTIME = 3;
    time = 0:Ts:RUNTIME;
    t_span = length(time);
    % traj. param for static sim
    waypoints = [0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0]';
    t_tp = [0,1,2];
    v_bound = zeros(6,3);
    fext = zeros(6,10); 
elseif sim_to_run == "ur5_dynamic_tauctrl.slx"
    % traj. for tau ctrl
    bspline = [1,4];
    % joint modes
    mode_unlock = zeros(1,6);
    mode_lock = ones(1,6);
    RUNTIME = 4;
    time= 0:Ts:4;
    t_span = length(time);
elseif sim_to_run == "ur5_dynamic_qctrl.slx"
    % shifting
    START = time(1);
    time_extended = 0:Ts:START;
    time = time-ones(1,t_span)*START;
    RUNTIME = time(t_span);
    % t_span = legnth(time)
end

% ext torque
TIM_MAG = 1/Ts;
delay = 0.5*TIM_MAG;
imp = RUNTIME-0.5;
tau_ext = 5;
fext = zeros(6,10); 
fext(:,10) = [0,0,0,0,0,0];


