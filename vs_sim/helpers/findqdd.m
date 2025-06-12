function [qdd] = findqdd(robot,u,q)

qin = q(1:6);
qdin = q(7:12);
textin = q(19:end);

qdd = (inv(massMatrix(robot,qin))*(u-gravityTorque(robot,qin)-velocityProduct(robot,qin,qdin)-textin));
end