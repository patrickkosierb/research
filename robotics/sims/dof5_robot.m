%% import robot from simulink %%
clc; close all;
[dof5_arm, arm_info] = importrobot('dof5_robotv2_block');

% additional details %
showdetails(dof5_arm)
% show(dof5_arm)

%% Export to URDF file
% exporter = urdfExporter(dof5_arm)
% writefile(exporter,OutputfileName="dof5_robot.urdf")

