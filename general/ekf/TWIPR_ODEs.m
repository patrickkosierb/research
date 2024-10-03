function y_dot = TWIPR_ODEs(t,y,T_motor,m_w,I_w,r, T_ext, m_b, I_b, l, g)
%INVERTED_PENDULUM_CART_ODES Summary of this function goes here
%   y = [ theta_L , theta_L_dot , theta_p , theta_p_dot ]' 
c_P = cos(y(3));
s_P = sin(y(3));

y_dot = zeros(4,1);
y_dot(1) = y(2);
y_dot(2) = ( r*c_P/( l+I_b/(l*m_b) ) * ( T_ext(t) + m_b*g*l*s_P) - l*r*m_b*y(4)^2*s_P + T_motor(t) ) / ( m_b*r^2 + 2*m_w*r^2 + 2*I_w - m_b*r^2*c_P^2/(1 + I_b/(l^2*m_b)) );
y_dot(3) = y(4);
y_dot(4) = ( T_ext(t) + l*r*m_b*y_dot(2)*c_P + m_b*g*l*s_P ) / (m_b*l^2 + I_b);
end