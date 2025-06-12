function [qdd] = findqdd_rtb(robot,u,q)

qin = [q(1),q(4),q(7),q(10),q(13),q(16)]';
qdin = [q(2),q(5),q(8),q(11),q(14),q(17)]';
textin = [q(3),q(6),q(9),q(12),q(15),q(18)]';

M = robot.inertia(qin');
G = robot.gravload(qin');
C = robot.coriolis(qin',qdin');
qdd = (inv(M)*(u-G'-C*qdin-textin));

end