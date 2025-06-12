clc;clear; close all;

% init

% index = 39418; 
index = 1;
% 31917 0.4V
% 39418 0.5V

% pre-processing
data = readtable('data/v250000_a50_121mv2_clean.csv'); 
t = data(index:end,1);
t = t{:,:};
x_meas = data(index:end,2);
x_meas = x_meas{:,:};
xdot = data(index:end,3);
xdot = xdot{:,:};
voltage = data(index:end,4);
voltage = voltage{:,:};
f_measured = data(index:end,5);
f_measured= f_measured{:,:};

% unit conversion to cm
x_meas = x_meas*100;
xdot = xdot*100;

% filtering noisey velocity
order = 3;
framelen = 11;
xdot = sgolayfilt(xdot,order,framelen);

state_est_bw
param_id_bw

trial = 3;

% Define the CSV file name
csvFilename = 'results.csv';

% Create the table
results = table(trial, SE_RT, PID_RT, rmse_x_se, rmse_xdot_se, rmse_x_pid, rmse_xdot_pid);

% Check if the file exists
if isfile(csvFilename)
    % Read the existing data
    oldData = readtable(csvFilename);

    % Append new results
    updatedData = [oldData; results];

    % Write the updated data back to CSV
    writetable(updatedData, csvFilename);
else
    % File does not exist, create a new CSV with headers
    writetable(results, csvFilename);
end

disp("Results saved and appended successfully!");