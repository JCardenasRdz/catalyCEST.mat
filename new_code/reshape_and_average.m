function [ Images3D ] = reshape_and_average( Images4D, cest_saturation_array )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% reshape and average     
[rows,cols, num_images ] = size(Images4D);  
num_offsets = length(cest_saturation_array);

Images4D = reshape(Images4D, ...
                     [rows,cols,num_offsets,  num_images / num_offsets  ]);        

Images3D  = mean( Images4D, 4); 

end

