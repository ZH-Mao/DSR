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

function [J_x_pose, J_x_intensity, I_Projec] = Fun_Jacobian(I0_xyz, I_grid,...
                                               I_x, I_y, I_z, i, T, I_mask, NumOfFrames)

R = T(1:3,1:3);
t = T(1:3,4);
[L, W, H]  = size(I_mask);
mask_grid  = griddedInterpolant(I_mask,'linear','none');
Projec_Local = I0_xyz*R' + kron(ones(size(I0_xyz,1),1),t');

index1 = find(Projec_Local(:,1)>=1 & Projec_Local(:,1)<=L &...
              Projec_Local(:,2)>=1 & Projec_Local(:,2)<=W &...
              Projec_Local(:,3)>=1 & Projec_Local(:,3)<=H );
Projec_Local = Projec_Local(index1,:);
isWhite      = mask_grid(Projec_Local);
index2       = find(isWhite==255);
Projec_Local = Projec_Local(index2,:);
index_j      = index1(index2);
num          = length(index2);

% I0(index_j) = 1;
% J_x_intensity = (I0(:))';
J_x_intensity = sparse(1:num,index_j,ones(num,1), num, length(I0_xyz));

a1 = Projec_Local(:,1); b1=Projec_Local(:,2);  c1=Projec_Local(:,3);
a2 = I_x(Projec_Local); b2=I_y(Projec_Local); c2=I_z(Projec_Local);
dPlocal_dPose = -[a2,b2,c2,b1.*c2-c1.*b2,c1.*a2-a1.*c2,a1.*b2-b1.*a2];

index_i    = zeros(1, NumOfFrames);
index_i(i) = 1;
J_x_pose   = kron(index_i, dPlocal_dPose);
J_x_pose   = sparse(J_x_pose);

I_Projec    = I_grid(Projec_Local);
end