clc;
close all;
clear

data = readtable(['data_voltagevarying_0volt.csv']);
T = 0.0147;
t = data(:,2);
t = t{:,:};
f_measured = data(:,6);
f_measured= f_measured{:,:};

x0 = data(:,3);
x0 = x0{:,:};
x1 = data(:,4);
x1 = x1{:,:};

f_lls = zeros(1,length(t));

c= 1155.0004257718545;
k= 308.6273993818113;
beta= 658.3020536059036;
gamma= -35.57841706427401;

z = 0;
n = 2;
for k = 1:length(t)-1
    dzdt = (x1(k)-gamma*abs(x1(k))*z*abs(z)^(n-1)-beta*x1(k)*abs(z)^n);
    f_lls(k) = c*x1(k)+k*z;
    z = z+dzdt*T;
end

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