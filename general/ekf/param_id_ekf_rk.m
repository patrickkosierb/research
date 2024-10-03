clc;clear;close all;
% system parameters
n = 5;
m = 2;
g = 9.81;
T = 0.1;
tspan = 0:T:10;

% sb-robot characteristics
m_b = 1;
h = 1.5;  
w = 0.4; 
d = 0.4; 
I_b = m_b * ( (d^2+h^2)/12 + d^2/4 )
l = h/2;
m_w = 0.3;
r = 0.3;
I_w = m_w/2*r^2; %Disk

% inputs
T_ext_input = @(t) 0;
T_motor_input = @(t) 0.5*sin(t); %* (heaviside(t - 1.2) - heaviside(t - 2.5));

% initial conditions
theta_L_0 = 0;
theta_L_dot_0 = 0;
theta_P_0 = 0;
theta_P_dot_0 = 0;
x0 = [ theta_L_0, theta_L_dot_0 , theta_P_0 , theta_P_dot_0]';

% simulating true dynamics
[t_actual,x_actual] = ode45(@(t,y) TWIPR_ODEs(t,y,T_motor_input,m_w,I_w,r, T_ext_input, m_b, I_b, l, g) , tspan , x0);
n_iter = length(t_actual);

% TRUE INERTIA
I_b_true = ones(n_iter,1)*I_b;

syms t T_motor theta_L theta_dot_L theta_P theta_dot_P I_b

X = [theta_L;theta_dot_L;theta_P;theta_dot_P;I_b];
u = [T_motor];

c_P = cos(theta_P);
s_P = sin(theta_P);

% kalman step up (predicting I_b)
Q = 10^-2*eye(n);
R = 10^-1*eye(m);

% State Transtion Function (Non-linear)
theta_L_double_dot_equation = ( r*c_P/( l+I_b/(l*m_b) ) * (m_b*g*l*s_P) - l*r*m_b*theta_dot_P^2*s_P + T_motor ) / ( m_b*r^2 + 2*m_w*r^2 + 2*I_w - m_b*r^2*c_P^2/(1 + I_b/(l^2*m_b)) );
theta_P_double_dot_equation = ( l*r*m_b*theta_L_double_dot_equation*c_P + m_b*g*l*s_P ) / (m_b*l^2 + I_b);

f = [theta_L; 
     theta_dot_L;
     theta_P;
     theta_dot_P;
     0]; % The parameter is constant

h = [theta_P+theta_dot_P*T,theta_dot_P+theta_P_double_dot_equation*T];

% jacobians
F = jacobian(f,X);
F_disc = eye(size(F))+F*T;
H = jacobian(h,X);
H_disc = eye(size(H))+H*T;

% making noisy data
measurement_noise = [mvnrnd(zeros(1,m),R,n_iter)];
y_noisy = [x_actual(:,3),x_actual(:,4)] + measurement_noise;

% Initial state estimates
x0 = [0; 0; 0; 0; 1]; % Initial state estimate

x_est = zeros(n, n_iter);
x_est(:,1) = x0;
x_pred = x_est;
P = Q;

% Simulation loop
for k = 1:n_iter

    u_input = [T_motor_input(t_actual(k))];
        
    % Prediction step
    % x_pred = vpa(subs(f,[X;u],[x_est(:,k);u_input]),3); % Predict the next state    

    % runge-kutta method
    k_1 = vpa(subs(f,[X;u],[x_est(:,k);u_input]),3);
    k_2 = vpa(subs(f,[X;u],[x_est(:,k)+0.5*T*k_1;T_motor_input(t_actual(k)+0.5*T)]),3);  
    k_3 = vpa(subs(f,[X;u],[x_est(:,k)+0.5*T*k_2;T_motor_input(t_actual(k)+0.5*T)]),3);  
    k_4 = vpa(subs(f,[X;u],[x_est(:,k)+T*k_3;T_motor_input(t_actual(k)+T)]),3);  
    x_pred = x_est(:,k)+(1/6)*(k_1+2*k_2+2*k_3+k_4)*T;  
   
    F_eval = vpa(subs(F_disc,[X;u],[x_est(:,k);u_input]),3);
    P_pred = F_eval*P*F_eval' + Q;

    % Update step
    H_eval = vpa(subs(H_disc,[X;u],[x_est(:,k);u_input]),3); % Evaluate Jacobian at the predicted state
    K = vpa(P_pred*H_eval'/(H_eval*P_pred*H_eval' + R),3); % Kalman gain    
    x_est(:,k+1) = x_pred + K*(y_noisy(k,:) - vpa(subs(h,[X(3:5);u],[x_est(3:5,k);u_input]),3))'; % Updated state estimate
    P = (eye(n) - K*H_eval)*P_pred; % Updated estimate covariance

end

I_b_avg = ones(n_iter,1)*mean(x_est(5,:));

% Plot the results
figure;
plot(tspan, x_est(5,1:n_iter));
hold on;
plot(tspan, I_b_true(1:n_iter)',"red");
hold on;
plot(tspan, I_b_avg(1:n_iter)',"green");
legend("I_b estimate", "I_b actual", "I-b avg est")

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
legend("Theta p estimate", "Theta p actual")



