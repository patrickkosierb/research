function dxdt = f_ODE_bw(t,x)

    beta = 130.79645522568575;
    gamma = -0.89773144490448;
    n_pa = 2;
    k = -60.078547581405964;
    c= 1201.3445613819722;

    dxdt = [
        x(2);
        c*x(2)+k*x(3);
        x(2)-beta*abs(x(2))*power(abs(x(3)),n_pa-1)*x(3)-gamma*x(2)*abs(x(3))^n_pa
        ];
    end
    