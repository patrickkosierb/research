%% ur5_sim with kf %%
init
run_model

%% kalman filter %%

% params
jnts = 1;
n=4;
m=jnts;

Q = 10e-3*eye(n*jnts); %e-6
R = 10e-9*eye(m); %e-4
P = Q;

% add noise to position 
q_noisey = noisey(q_true,m,R,t_len);

% init est
x_est = zeros(n*jnts, t_len);

m = 1;
cor = m;
g = 0;
A = [ 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; %q
      0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; %qd
      0 m*cor 0 m 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; %qdd
      0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
];
B = [0 0 1 0];
C = [1 0 0 0];

% kf loop
tic;
for k=1:t_len-1
% m = 1/(massMatrix(ur5e,x_est(1,k)));
% cor = 1;
% u = m*tau_jint_out(:,k)-m*tau_G_out(:,k);
% [x_est(:,k+1), P] = kf(x_est(:,k), q_noisey(1,k), u, P, A, B, C, Q, R, n*jnts) ;
end
toc
