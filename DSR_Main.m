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

clear;
clc;
close all;
format long;

addpath(genpath('./dataset/'))
filename   = dir('./dataset/*.mat');
NumOfFrames= size(filename,1);

I          = cell(NumOfFrames,1); % Local 3D images
mask       = cell(NumOfFrames,1); % Mask matrices to mask the invalide voxels
I_grid     = cell(NumOfFrames,1); % Trilinear interpolation of local images
I_gradient = cell(NumOfFrames,3); % Gradient of local images
T          = cell(NumOfFrames,1); % Global Poses of local images

% Initial guesses of poses, obtained from the pairwise registration
Li         = importdata('RelativePoses_Pairwise.mat'); 

for i = 1:NumOfFrames
    I{i}            = double(importdata(filename(i).name));
    mask{i}         = Fun_mask(I{i}, 10);
    I_grid{i}       = griddedInterpolant(I{i},'linear','none');
    [Ix,Iy,Iz]      = gradient(I{i});
    I_gradient{i,1} = griddedInterpolant(Ix,'linear','none');
    I_gradient{i,2} = griddedInterpolant(Iy,'linear','none');
    I_gradient{i,3} = griddedInterpolant(Iz,'linear','none');
    if i ==1
        T_tmp = eye(4);
        T{i}  = T_tmp;        
    else
        T_tmp = Lie2T(Li(:,i-1))*T_tmp;
        T{i}  = T_tmp;
    end
end

I0 = zeros(size(I{1}));  % Create a grid of a panoramic image

clear I Ix Iy Iz T_tmp

[L0,W0,H0] = size(I0);
[x0,y0,z0] = ind2sub([L0,W0,H0],(1:L0*W0*H0)');
I0_xyz     = [x0,y0,z0]; % Coordinates of voxels in the panoramic image

step         = 1;  % Initialize step length of poses
m            = 0;  % Initialize iteration counter
Iterations   = []; % Save the iteration number
D_length     = []; % Save step length of poses in each iteration
Num_voxels   = []; % Save number of intensity differences in each iteration
time         = []; % Computational cost for every iteration
maxIteration = 200;

while step >1e-10 && m < maxIteration
    tic;
    m = m+1;
    Iterations = [Iterations,m];
    
    J_x_Poses       = [];
    J_x_Intensities = [];
    I_Projec        = [];
    for i = 1:NumOfFrames
        % Calculate Jacobian matrix
        [J_x_pose, J_x_intensity, I_Projec_tmp] = Fun_Jacobian(I0_xyz, I_grid{i},...
            I_gradient{i,2}, I_gradient{i,1}, I_gradient{i,3}, i, T{i}, mask{i},NumOfFrames);
        J_x_pose(:,1:6) = [];
        J_x_Poses       = [J_x_Poses;J_x_pose];
        J_x_Intensities = [J_x_Intensities;J_x_intensity];
        I_Projec        = [I_Projec;I_Projec_tmp];
    end
    
    % Components of Hessian matrix in GN equation
    H_PP  = (J_x_Poses)'*(J_x_Poses);
    H_PI  = (J_x_Poses)'*(J_x_Intensities);
    H_II  = (J_x_Intensities)'*J_x_Intensities;
    
    I_Projec  = sparse(I_Projec);
    Num_voxels=[Num_voxels, size(I_Projec,1)];
    [i,j]     = find(H_II);
    H_II_dig  = nonzeros(H_II);
    inv_H_II  = sparse(i,j,1./H_II_dig,size(H_II,1),size(H_II,2));
    
    % Step change of poses
    d_pose    = -(H_PP-H_PI*inv_H_II*(H_PI)')\(H_PI*inv_H_II*J_x_Intensities'*I_Projec-J_x_Poses'*I_Projec);

    % Update poses
    for k = 1: NumOfFrames-1
        dT     = Lie2T(d_pose(1+6*(k-1):6*k));
        T{k+1} = dT*T{k+1};
    end
    step     = norm(d_pose);
    D_length = [D_length,step];
    time     = [time,toc];
end
% Plot the iteration process to confirm if the code converges
plot(Iterations,D_length,'ro-');
grid on;
xlabel('Iteration');
ylabel('Step Length');
title('The Step Changes of Poses');

save Poses_DSR.mat T Iterations D_length time Num_voxels