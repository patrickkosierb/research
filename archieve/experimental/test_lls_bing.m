clc;
close all;
clear

data = readtable(['data_voltagevarying_0volt.csv']);

t = data(:,2);
t = t{:,:};
f_measured = data(:,6);
f_measured= f_measured{:,:};

x0 = data(:,3);
x0 = x0{:,:};
x1 = data(:,4);
x1 = x1{:,:};

order = 8;
framelen = 31;

x1 = sgolayfilt(x1,order,framelen);

f_lls = zeros(1,length(t));
% 
% cd= 271.2174064452846;
% fc= 5.984673215294906;
% f0= -8.274894565102125;

% cd= 271.2174071412453;
% fc= 5.98467321759326;
% f0= -8.27489455281873;
cd= 366.0735256598746
f= 5.451537602810378
f0= -8.273570727280303
% Matlab
% cd = 352.5194;
% fc = 5.5396; 
% f0 = -8.2707;

for k = 1:length(t)-1
    f_lls(k) = cd*x1(k)+fc*sign(x1(k))+f0;
end

squaredDifferences = (f_measured' - f_lls).^2;
% Compute the mean of squared differences
mse = mean(squaredDifferences);
% Compute the RMSE by taking the square root of MSE
rmse = sqrt(mse)

% rmse_force = sqrt(mean(f_measured' - f_lls,2).^2)% this is just vectorized

figure;
plot(t, f_measured); 
hold on 
plot(t,f_lls);
title("Force Vs. Time");
legend("Measured","LLS");
xlabel('Time(s)')
ylabel('F')
hold off