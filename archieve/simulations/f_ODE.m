function dxdt = f_ODE(t,x)

factor = 2;

omega = 0.5;
c = 0.3*factor;
k = 1;
beta = 2*factor;
gamma = 1;
n_pa = 2;

dxdt = [
    x(2);
    5*sin(omega*t)-(c*x(2)+k*x(3));
    x(2)-beta*abs(x(2))*power(abs(x(3)),n_pa-1)*x(3)-gamma*x(2)*abs(x(3))^n_pa
    ];
end
