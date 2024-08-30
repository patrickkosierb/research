function dxdt = boucwen(t,x)

    c0= 11.825884573691747;
    k0= 0.6168415593308919;
    alpha= 44.45425042190868;
    beta= -3061.832988105183;
    gamma= 3078.384741848536;
    A= 2.761665766402648;
    x0= 7.191401303755131;
    n_po = 2;
    period = 15;
    
    w = 2*pi*1/period;

    dxdt = [
        x(2);
        100*0.0375*w^2*sin(w*t-pi/2)-(c0*x(2)+k0*(x(1)-x0)+alpha*x(3));
        -gamma*abs(x(2))*x(3)*abs(x(3))^(n_po-1)-beta*x(2)*abs(x(3))^n_po+A*x(2);
        ];
    end