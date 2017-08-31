function [ mask_tumor, mask_noise] = select_ROIS( imgs3D)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% Normalize
imgs3D = imgs3D ./ max(imgs3D(:));
%% SELECT ROIs
I  = max(imgs3D,[],3);
% tumor
imagesc(I); colormap('jet');
title(strcat('Choose TUMOR ROI'));
mask_tumor=roipoly; 

% noise
imagesc(I); colormap('parula');
title(strcat('Choose NOISE ROI'));
mask_noise=roipoly; close all;

end

