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
% The code is used to generate a mask matrix which main area is white(255)
% Image: input image
% m: exclude m-pixel-width on the image boundary
% ===============================================================================================


function I=Fun_mask(Image, m)

[L,W,H]=size(Image);
Image(1,:,:)=0;
Image(L,:,:)=0;
Image(:,1,:)=0;
Image(:,W,:)=0;
Image(:,:,1)=0;
Image(:,:,H)=0;
res=Image;
if m==0
    res(Image~=0)=255;
    I=double(res);
else
    for n=1:m
        for k=2:H-1
            for j=2:W-1
                for i=2:L-1
                    if Image(i,j,k)==0
                        continue;
                    end
                    if (Image(i-1,j,k)==0)||(Image(i+1,j,k)==0)||(Image(i,j-1,k)==0)||(Image(i,j+1,k)==0)||(Image(i,j,k-1)==0)||(Image(i,j,k+1)==0)
                        res(i,j,k)=0;
                    end
                end
            end
        end
        Image=res;
    end
    I=double(res);
    I(I~=0)=255;
end
end