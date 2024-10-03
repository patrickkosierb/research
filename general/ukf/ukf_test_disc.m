clc; clear; close all;
syms pos pos_dot theta theta_dot

% sim param
T = 0.001;
t_span = 0:T:10;
t_len = length(t_span);
n = 4;
m = 1;

g = 9.81;
l = 1;
mass = 10;
k=1;

% Initialize state vector
X = [pos;pos_dot;theta;theta_dot]
X_EQ = [ 0; 0; 0; 0];

f = [pos_dot;
    (l+pos)*theta_dot^2+g*cos(theta)-k/mass*pos;
    theta_dot;
    (-2*pos_dot*theta_dot-g*sin(theta))/(l+pos);
]

% Finding A, B, C (to linearize, we need to use the jacobian)
A = jacobian(f,X);
A_eq = subs(A,X,X_EQ);
C = [0 0 1 0]; % reading theta
D = [0];

A_dis= eye(n)+T*A_eq;
x_dot = zeros(n,t_len);
y = zeros(m,t_len);

% initial conditions
x_dot(:,1) = [0.1; 0; pi/2; 0];
for k=1:t_len-1
    x_dot(:,k+1) = A_dis*x_dot(:,k);
    y(k+1) = C*x_dot(:,k+1);
end

x0 = [0.1; 0; pi/2; 0];

%% UKF 
f_ukf = @(x,u)[
    x(1)+T*x(2);
    x(2)+T*((x(1)+l)*x(4)^2-k/m*x(1)+g*cos(x(3)))
    x(3)+T*x(4);
    x(4)+T*(-g*sin(x(3))-2*x(2)*x(4))/(x(1)+l);
];

h_ukf = @(x,u)[x(3)+T*x(4);];

Q = 1e-6*eye(n); 
R = 1e-5*eye(m); 
P_ukf = Q;

% making noisy data
theta_actual = y;
noise = [mvnrnd(zeros(1,m),R,t_len)]';
theta_noise = [theta_actual] + noise;

x_est = zeros(n, t_len);
x_est(:,1) = x0;
theta_est = zeros(1,t_len);

% ukf loop
for k = 1:t_len-1
    [x_est(:,k+1), P_ukf] = ukf(x_est(:,k), theta_noise(k), 0, P_ukf, f_ukf, h_ukf, Q, R);
    theta_est(k+1) = h_ukf(x_est(:,k+1),0);
end

rmse_x=sqrt(mean(x_est(1,:)-x_dot(1,:))^2)
rmse_theta=sqrt(mean(theta_est-theta_actual)^2)

figure;
plot(t_span,theta_actual)
hold on;
plot(t_span,theta_est)
legend("Actual","UKF"  )

figure;
plot(t_span,x_dot(1,:))
hold on;
plot(t_span,x_est(1,:))
legend("Actual ode","UKF")
