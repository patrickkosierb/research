function x_dot = pendODEs(t,x,m,l,k,g)

% for ease
c_P = cos(x(3));
s_P = sin(x(3));

% define state vector (this is also the state transition function)
x_dot = zeros(4,1);

x_dot(1) = x(2);
x_dot(2) = (l+x(1))*x(4)^2+g*c_P-k/m*x(1);
x_dot(3) = x(4);
x_dot(4) = (-2*x(2)*x(4)-g*s_P)/(l+x(1));

end