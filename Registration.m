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


function I2 =Registration(I1, T)
I1 = double(I1);

I1grid=griddedInterpolant(I1,'linear','none');
[l,w,h]=size(I1); 
I2 = zeros(l,w,h);
R  = T(1:3,1:3);
R_ = inv(R);
t  = T(1:3,4);

[x2,y2,z2] = ind2sub([l,w,h],(1:l*w*h)');
I2_xyz     = [x2,y2,z2];
P1         = (I2_xyz-kron(ones(l*w*h,1),t'))*R_';

index1 = find(P1(:,1)>=1 & P1(:,1)<=l &...
              P1(:,2)>=1 & P1(:,2)<=w &...
              P1(:,3)>=1 & P1(:,3)<=h );
P1     = P1(index1,:);
index2 = find(I1grid(P1));
P1     = P1(index2,:);
P2     = I2_xyz(index1(index2),:);
I2(index1(index2)) = I1grid(P1);
end
