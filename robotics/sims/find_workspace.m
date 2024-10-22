%% find robot workspace %%
clc; close all; clear;

% import robot into rigidbody
[dof5_arm, arm_info] = importrobot('dof5_robotv2_block');

% workspace param init
n_points = 25;
joint_lim = linspace(-pi/2, pi/2, n_points);
workspace = [];

% find all possible joint combinations (1 loop per joint)
for q = 1:n_points
    for n = 1:n_points
        for k = 1:n_points
            % Create a configuration vector for the robot
            config = homeConfiguration(dof5_arm);
            config(1).JointPosition = joint_lim(q);
            config(2).JointPosition = joint_lim(n);
            config(3).JointPosition = joint_lim(k);

            % Compute the forward kinematics for the current configuration
            endEffectorTform = getTransform(dof5_arm, config, 'Body6');
            position = tform2trvec(endEffectorTform);

            % Collect the end-effector position
            workspace = [workspace; position];

        end
    end
end

% plot the workspace
plot3(workspace(:,1), workspace(:,2), workspace(:,3));
% scatter3(workspace(:,1), workspace(:,2), workspace(:,3), '.');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Robot Workspace');
grid on;

