% Script to import sm_universalUR5.urdf
% Copyright 2021-2024 The MathWorks, Inc.

%% Import URDF file to create Simscape Multibody model
addpath([pwd filesep 'Geometry'])
[mdl_h] = smimport('sm_universalUR5.urdf','ModelName','ur5_block');
mdl_name = getfullname(mdl_h);

%% Update diagram, note initial robot position
set_param(mdl_h,'SimulationCommand','update')

