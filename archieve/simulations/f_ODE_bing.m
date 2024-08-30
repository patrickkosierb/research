function dxdt = f_ODE_bing(t,x)

omega = 0.5;
cd = 0.3;
fc = 1;

dxdt = [
    x(2);
    5*sin(omega*t)-(cd*x(2)+fc*sign(x(2)))
    ];
end
