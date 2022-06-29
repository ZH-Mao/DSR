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
% ===============================================================================================
% This code is used to estimate relative pose between Frame I and Frame I+1;
% The results are saved as Lie algebra;
% ===============================================================================================

clear;
close all;
clc;

format long;
addpath('./dataset')

filename= dir('./dataset/*.mat');
fn      = size(filename,1);
Li     = zeros(6,fn-1);

for u = 1:fn-1
    I1name = filename(u).name;
    I1 = double(importdata(I1name));
    I2name = filename(u+1).name;
    I2 = double(importdata(I2name));
    
    I1_grid=griddedInterpolant(I1,'linear','none');
    I2_grid=griddedInterpolant(I2,'linear','none');
    
    [I2x,I2y,I2z]=gradient(I2);
    I2_x=griddedInterpolant(I2x,'linear','none');
    I2_y=griddedInterpolant(I2y,'linear','none');
    I2_z=griddedInterpolant(I2z,'linear','none');
    [l_,w_,h_]=size(I2);
    
    mask1=Fun_mask(I1,10);
    mask2=Fun_mask(I2,10);
    
    mask2_grid=griddedInterpolant(mask1,'linear','none');
    
    index      = find(mask1~=0);
    [l,w,h]    = size(I1);
    [x1,y1,z1] = ind2sub([l,w,h],index);
    I1_xyz     = [x1,y1,z1];
    
    d         = ones(6,1);
    m         = 1;
    T12       = eye(4);
    maxIteration = 300;
    
    while (norm(d)>1e-10) && (m<=maxIteration)
        R = T12(1:3,1:3);
        t = T12(1:3,4);
        P2 = I1_xyz*R' + kron(ones(size(I1_xyz,1),1),t');
        
        index1 = find(P2(:,1)>=1 & P2(:,1)<=l_ &...
            P2(:,2)>=1 & P2(:,2)<=w_ &...
            P2(:,3)>=1 & P2(:,3)<=h_ );
        P2      = P2(index1,:);
        isWhite = mask2_grid(P2);
        index2  = find(isWhite==255);
        P2      = P2(index2,:);
        P1      = I1_xyz(index1(index2),:);
        num     = size(P1,1);
        e       = I1_grid(P1)- I2_grid(P2);
        a1 = P2(:,1); b1=P2(:,2);  c1=P2(:,3);
        a2 = I2_y(P2); b2=I2_x(P2); c2=I2_z(P2);
        J  = -[a2,b2,c2,b1.*c2-c1.*b2,c1.*a2-a1.*c2,a1.*b2-b1.*a2];
        e  = sparse(e);
        J  = sparse(J);
        d =(J'*J)\(-J'*e);
        dT = Lie2T(d);
        T12  = dT*T12;
        objvalues  = e'*e;
        steplength = norm(d);
        
        disp(m);
        disp('Values of objective function is:');
        disp(objvalues);
        disp('Step length is:');
        disp(steplength);
        m=m+1;
    end
    Li(:,u) = T2Lie(T12);
end
savefile = sprintf('RelativePoses_Pairwise.mat');
save( savefile,'Li'); 
