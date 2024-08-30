dt=0.05;
t = 0:dt:30;

T = 15;
w = 2*pi*1/T;
x = -0.0375*sin(w*t-pi/2)-0.0375;

xdot = -0.0375*w*cos(w*t-pi/2);
xddot = 0.0375*w^2*sin(w*t-pi/2);

plot(t,x);
hold all
plot(t,xdot);
plot(t,xddot)

