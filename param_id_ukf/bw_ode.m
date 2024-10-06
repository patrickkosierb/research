function x_dot = bw_ode(c0,k0,alpha,beta,gamma,A,x0)
x_dot = zeros(3,1);

z = 1; %evolutionary variable

x_dot = [x(2)k; %position
         1; %velocity
         c0*x(2)+k0*(x(1)-x0)+alpha*z; % force
        ];



end