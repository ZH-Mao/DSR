% ===============================================================================================
% DSR: Direct Simultaneous Registration for Multiple 3D Images
% Version: 1.0
% ===============================================================================================
% 
% Copyright (C) 2022 Zhehua Mao, Liang Zhao, and Shoudong Huang
% University of Technology, Sydney, Australia
% 
% Authors:  Zhehua Mao         -- Zhehua.Mao@student.uts.edu.au
%           Liang Zhao         -- Liang.Zhao@uts.edu.au 
%           Shoudong Huang     -- Shoudong.Huang@uts.edu.au
% 
%           Robotics Institute
%           Faculty of Engineering and Information Technology
%           University of Technology, Sydney
%           NSW 2007, Australia
% 
% Please contact Zhehua Mao {Zhehua.Mao@student.uts.edu.au} if you have any questions about the code.



function  x=T2Lie(T)
R        = T(1:3,1:3);
t        = T(1:3,4);
theta    = acos((trace(R)-1)/2);
temp_R_  = theta/(2*sin(theta))*(R-R');
R_       = [temp_R_(3,2); temp_R_(1,3); temp_R_(2,1)];
theta    = norm(R_);
n        = R_/theta;
n_       = [0 -n(3) n(2);n(3) 0 -n(1);-n(2) n(1) 0];
J        = (sin(theta)/theta)*eye(3)+(1-sin(theta)/theta)*(n*n')+((1-cos(theta))/theta)*n_;
t_       = J\t; 
x        = [t_;R_];
end