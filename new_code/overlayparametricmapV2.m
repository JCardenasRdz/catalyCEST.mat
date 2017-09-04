function [] = overlayparametricmapV2( AnatomicalImage,MAP,mask ,range, titulo)
%UNTITLED Summary of this function goes here
%   [] = overlayparametricmap( AnatomicalImage,MAP,mask )
% AnatomicalImage = grey anatonical reference
% MAP= parametric map
% mask = binary mask

% GFC version 1.0
% Sep, 2017 !!

% parametric
mymask = jet(64);
%mymask(1,:) = 0;
%mymask(1,:)=0;


% check sizes
[rows,cols] = size(AnatomicalImage);

% resize parametric map
MAP = imresize(MAP,[rows cols]);

% resize mask
mask = imresize(mask,[rows cols]);

baseImage=AnatomicalImage;
parameterROIImage=MAP;

baseImage = baseImage/(max(baseImage(:))); % normalize base (anatomical) image
rgbSlice  = baseImage(:,:,[1 1 1]);        % converting to RGB (ignore colormaps)
figure();
           % show parametric image

    imshow(parameterROIImage, range);  

colormap(mymask);                           % apply colormap
hold on;
h = imshow(rgbSlice);                      % superimpose anatomical image
set(h, 'AlphaData', mask==0);      % make pixels in the ROI transparent
colorbar; 
title(titulo)
end

