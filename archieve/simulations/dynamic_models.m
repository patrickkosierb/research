clear;
clc;
close all;

T=0.05;
t = 0:T:30;

omega = 1;
x = 5*sin(omega*t);
xdot = 5*omega*cos(omega*t);

c0 = 300;
fc = 30;
f0 = -5;

fmr = zeros(1,length(t));

for k = 1:length(t)-1
    fmr(k+1) = c0*xdot(k)+fc*sign(xdot(k))+f0;
end

plot(t,x)
hold on 
plot(t,xdot)
hold off

figure;
tiledlayout(2,1)
nexttile
plot(t,fmr)
nexttile
plot(fmr,xdot)

