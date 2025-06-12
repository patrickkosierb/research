%% Kalman filter
function [x, P] = kf(x, z, u, P, A, B, C, Q, R, n)
x = A*x + B*u; 
P = A*P*A' + Q;
K = P*C'/(C*P*C' + R); 
x = x + K*(z - C*x);
P = (eye(n) - K*C)*P*(eye(n) - K*C)' + K*R*K';
end