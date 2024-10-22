%% run find_workspace.m if robots workspace isn't defined yet %%
clc; close all; 

% test trajectory from robot workspace
t_span = linspace(0,10,25);
traj1 = [workspace(6900:6924,1)'; workspace(6900:6924,2)'; workspace(6900:6924,3)'];
traj1_v = zeros(3, 25);

% forward kinenmatics output from simulation
fwd_x = out.fwd_out.data(1,:);
fwd_y = out.fwd_out.data(2,:);
fwd_z = out.fwd_out.data(3,:);
fwd = [fwd_x; fwd_y; fwd_z];

figure;
plot3(traj1(1,:)', traj1(2,:)', traj1(3,:)'); hold on;
plot3(fwd_x', fwd_y', fwd_z');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Test Trajectory');
legend("Desired", "Actual")
grid on;

% % random waypoints
% wp1 = [2.5,2.5,2.5;1.75,1.75,1.75;0,0,0];
% wp2 = [1,1,1;1,1,1;1,1,1];
% wp3 = [2,2,2;2.5,2.5,2.5;0,0,0];
% wp4 = [2,3,2.5;1,1,1;0,1,2];

% size(traj1)
% size(fwd)
% RMSE = [];
% for i=1:3
% 
%     comp = sqrt(mean((traj1-fwd).^2))
%     RMSE = [RMSE; comp]
% end
% RMSE
