
% inspired by https://github.com/petercorke/robotics-toolbox-matlab/blob/master/models/mdl_ur5.m#L55
function r = mdl_ur5e()
    
    deg = pi/180;
    
    % robot length values (metres)
    a = [0, -0.42500, -0.39225, 0, 0, 0]';

    d = [0.1625, 0, 0, 0.1333, 0.0997, 0.0996]';

    alpha = [1.570796327, 0, 0, 1.570796327, -1.570796327, 0]';
    
    theta = zeros(6,1);
    
    DH = [theta d a alpha];

    mass = [3.7000, 8.3930, 2.33, 1.2190, 1.2190, 0.1897];

    center_of_mass = [
        0,-0.02561, 0.00193
        0.2125, 0, 0.11336
        0.15, 0, 0.0265
        0, -0.0018, 0.01634
        0, -0.0018, 0.01634
        0, 0, -0.001159];
    
    % and build a serial link manipulator
    
    % offsets from the table on page 4, "Mico" angles are the passed joint
    % angles.  "DH Algo" are the result after adding the joint angle offset.

    robot = SerialLink(DH, ...
        'name', 'UR5e', 'manufacturer', 'Universal Robotics');
    
    % add the mass data, no inertia available
    links = robot.links;
    for i=1:6
        links(i).m = mass(i);
        links(i).r = center_of_mass(i,:);
        links(i).I = [0,0,0,0,0,0];
        links(i).Jm = 0;
    end
    
    % links(1).I = [0.010267495893, 0.010267495893, 0.00666, 0,0,0];
    % links(2).I = [0.22689067591, 0.22689067591, 0.0151074, 0,0,0];
    % links(3).I = [0.049443313556, 0.049443313556, 0.004095, 0,0,0];
    % links(4).I = [0.111172755531, 0.111172755531, 0.21942, 0,0,0];
    % links(5).I = [0.111172755531, 0.111172755531, 0.21942, 0,0,0];
    % links(6).I = [0.0171364731454, 0.0171364731454, 0.033822, 0,0,0];

    % place the variables into the global workspace
    if nargin == 1
        r = robot;
    elseif nargin == 0
        assignin('caller', 'ur5e', robot);
        assignin('caller', 'qz', [0 0 0 0 0 0]); % zero angles
        assignin('caller', 'qr', [180 0 0 0 90 0]*deg); % vertical pose as per Fig 2
    end
end