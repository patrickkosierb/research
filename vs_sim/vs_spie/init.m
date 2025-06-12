%% init virtual sensor sim %%
clear; clc; close all;

addpath(pwd+"/helpers/")
addpath(pwd+"/Geometry/")

% simulation param 
ur5_rbt = importrobot("ur5_block",DataFormat="column");
mdl_ur5
runtime = 3;
Ts = 0.01;
t = 0:Ts:runtime;
t_len = length(t);
TIM_MAG = 1/Ts;
% ext torque 
delay = 0.5*TIM_MAG;
imp = 2.5;
tau_ext = 5;  %0.1 or 5
% traj. param 
waypoints = [0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0]';
t_tp = [0,1,2];
v_bound = zeros(6,runtime);
fext = zeros(6,10); 




