run_model
run_ukf
run_fomo
run_gm_kf

figure;
plot(t,tau_ext_ideal(6,:), 'b');
hold on;
plot(t,x_est(24,:),'r')
plot(t,-r(6,:),'g');
plot(t,-x_estgm(12,:),'m');
xlabel("Time(s)")
ylabel("\tau_{ext}(Nm)")
legend("Actual", "FOMO", "GMKF")
title("Tool Z-axis \tau_{ext}")


trial = "-1";

% Define file names
figFilename = 'vs_spie/results_wgmkf/trial'+trial+'.fig';  % MATLAB Figure File
pngFilename = 'vs_spie/results_wgmkf/trial'+trial+'.png';  % PNG Image File

% Save the figure
savefig(figFilename);   % Saves in MATLAB's .fig format
saveas(gcf, pngFilename);  % Saves as PNG image

% Define the CSV file name
csvFilename = 'vs_spie/results_wgmkf/testing.csv';

% Create the table
results = table(trial, UKF_RT, FOMO_RT, GMKF_RT, rmsetau_ukf, rmsetau_fomo, rmsetau_gmkf);

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