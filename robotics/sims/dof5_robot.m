clc; close all; clear;
[dof5_arm, arm_info] = importrobot('dof5_robot_block')


showdetails(dof5_arm)
show(dof5_arm)
% Export the URDF file
exporter = urdfExporter(dof5_arm)
writefile(exporter,OutputfileName="dof5_robot.urdf")

