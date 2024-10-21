%% ROBOT TEST BED (RTB) PARAMETERS 
clc;clear;

% setup enclosure (m) based on table measurements 
h = 1.8;
w = 0.70;
l = 1.15;
h_base = 0.2; % assuming some power or computer module under base frame
d_base = 0.2;

%canadarm3 specs (metric) from https://www.asc-csa.gc.ca/eng/iss/canadarm2/canadarm-canadarm2-canadarm3-comparative-table.asp
l_c3 = 8.5;
m_c3 = 715;
d_c3 = 0.23;

% calculating reachability for max length of rtb (i.e. length of 2 booms)
p_reachable = [w/2;l/2;h];
rtb_len_max = sqrt(sum(p_reachable.^2));
rtb_len_max_base = sqrt(p_reachable(1)^2+p_reachable(2)^2+(p_reachable(3)-h_base)^2);  %including base
rtb_len_max_base = rtb_len_max_base - 2/100; % 2cm margin

R_c3_to_rtb = rtb_len_max_base/l_c3 %scale for test bed approx 1/5th

% rtb parameters
l_rtb = rtb_len_max_base
m_rtb = R_c3_to_rtb*m_c3
d_rtb = R_c3_to_rtb*d_c3

d_rtb_joint = d_rtb+0.05 %made it up
l_rtb_ee = 0.15 %made it up
l_rtb_base = l_rtb_ee;

l_rtb_boom = (rtb_len_max_base-l_rtb_ee-l_rtb_base)/2


% Define the data for the table
% param = {'ICE RTB'; 'Jane'; 'Mary'; 'Tom'};
% ages = [28; 34; 22; 45];
% heights = [5.9; 5.6; 5.7; 6.0];  % in feet
% weights = [180; 135; 150; 200];   % in lbs
% 
% % Create the table
% T = table(names, ages, heights, weights);
% 
% % Add column names
% T.Properties.VariableNames = {'Name', 'Age', 'Height_ft', 'Weight_lbs'};
% 
% % Display the table
% disp(T);