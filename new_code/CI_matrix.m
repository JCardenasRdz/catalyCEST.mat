function [ CI_map ] = CI_matrix( CESTfit, mask_toi)
%UNTITLED Summary of this function goes here
% converts output from CEST MRI function into a paramtric map of the lower
% CI of a particular peak.
% Gaby uses this filter voxels
%   Detailed explanation goes here

num_of_voxels = length(CESTfit);
lower_ci_SA_peak = zeros(num_of_voxels,1);

for j = 1:num_of_voxels
    lower_ci_SA_peak(j)= CESTfit{j}.cfit.ci(2,1); 
    
end

% to create parametric map
[rows,cols] = size(mask_toi);
indices = find(mask_toi);

CI_map = nan(rows,cols);
CI_map(indices) = lower_ci_SA_peak;

end

