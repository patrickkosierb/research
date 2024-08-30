function dxdt = f_ODE_bing(t,x)

cd = 0.3;
fc = 1;

dxdt = [
    x(2);
    cd*x(2)+fc*sign(x(2))
    ];
end
