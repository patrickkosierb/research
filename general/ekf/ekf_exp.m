clc;clear;close all;

addpath("C:\Users\pk\Documents\ta\4ss3\data\");

% pre-processing
data = readtable('imu_data_4.txt');

t = data(:,1);
tspan = t{:,:}./1000-1.968;

theta_p_meas = data(:,2);
y_noisy = theta_p_meas{:,:}.*(pi/180);

% system parameters
n = 4;
m = 1;
g = 9.81;
T = 0.01;
% tspan = 0:T:2;

% sb-robot characteristics
m_b = 1.5;
h = 0.14;  
w = 0.14; 
d = 0.08; 
I_b = m_b * ( (d^2+h^2)/12 + d^2/4 );
l = h/4;
m_w = 0.1;
r = 0.035;
I_w = m_w/2*r^2; %Disk

% inputs
T_ext_input = @(t)  0;
T_motor_input = @(t) 0;

% initial conditions
theta_L_0 = 0;
theta_L_dot_0 = 0;
theta_P_0 = pi/6;
theta_P_dot_0 = 0;
x0 = [ theta_L_0, theta_L_dot_0 , theta_P_0 , theta_P_dot_0]';

% simulating true dynamics
[t_actual,x_actual] = ode45(@(t,y) TWIPR_ODEs(t,y,T_motor_input,m_w,I_w,r,T_ext_input, m_b, I_b, l, g) , tspan , x0);
n_iter = length(t_actual);

for k=1:n_iter
    if(x_actual(k,3)>pi/2)
        x_actual(k,3)=pi/2;
    end
end

syms t T_motor theta_L theta_dot_L theta_P theta_dot_P

X = [theta_L;theta_dot_L;theta_P;theta_dot_P];
u = [T_motor];

c_P = cos(theta_P);
s_P = sin(theta_P);

% kalman step
Q = 10^-4*eye(n);
R = 10^-3*eye(m);

% State Transtion Function (Non-linear)
theta_L_double_dot_equation = ( r*c_P/( l+I_b/(l*m_b) ) * (m_b*g*l*s_P) - l*r*m_b*theta_dot_P^2*s_P + T_motor ) / ( m_b*r^2 + 2*m_w*r^2 + 2*I_w - m_b*r^2*c_P^2/(1 + I_b/(l^2*m_b)) );
theta_P_double_dot_equation = ( l*r*m_b*theta_L_double_dot_equation*c_P + m_b*g*l*s_P ) / (m_b*l^2 + I_b);


f = [theta_L+theta_dot_L*T; 
     theta_dot_L+theta_L_double_dot_equation*T;
     theta_P+theta_dot_P*T;
     theta_dot_P+theta_P_double_dot_equation*T;
    ];

h = [theta_P+theta_dot_P*T];

% jacobians
F = jacobian(f,X);
H = jacobian(h,X);

% making noisy data

% measurement_noise = [mvnrnd(zeros(1,m),R,n_iter)];
% y_noisy = [x_actual(:,3)] + measurement_noise;

x_est = zeros(n, n_iter);
x_est(:,1) = x0;
x_pred = x_est;
P = Q;

% Simulation loop
for k = 1:n_iter

    u_input = [T_motor_input(t_actual(k))];
        
    % Prediction step
    x_pred = vpa(subs(f,[X;u],[x_est(:,k);u_input]),3); % Predict the next state
    F_eval = vpa(subs(F,[X;u],[x_est(:,k);u_input]),3);
    P_pred = F_eval*P*F_eval' + Q;
    
    % Update step
    H_eval = vpa(subs(H,[X;u],[x_est(:,k);u_input]),3); % Evaluate Jacobian at the predicted state
    K = vpa(P_pred*H_eval'/(H_eval*P_pred*H_eval' + R),3); % Kalman gain    
    x_est(:,k+1) = x_pred + K*(y_noisy(k,:) - vpa(subs(h,[X(3:4);u],[x_est(3:4,k);u_input]),3))'; % Updated state estimate
    P = (eye(n) - K*H_eval)*P_pred; % Updated estimate covariance


    if(x_est(3,k)>pi/2)
        x_est(3,k)=pi/2;
    end

end


figure;
tiledlayout(2,1);
nexttile;
plot(tspan, x_est(1,1:n_iter));
hold on;
plot(tspan, x_actual(1:n_iter,1)');
legend("Theta L estimate", "Theta L actual")
nexttile;
plot(tspan, x_est(3,1:n_iter));
hold on;
plot(tspan, x_actual(1:n_iter,3)');
plot(tspan, y_noisy(1:n_iter));
legend("Theta p estimate", "Theta p simulated", "Theta p measured")



