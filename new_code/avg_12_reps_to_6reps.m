function [ Data06reps ] = avg_12_reps_to_6reps( Data12reps )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[rows,cols,~] = size(Data12reps);
Data06reps = nan(rows,cols,6);


index_to_start = [1,3,5,7,9,11];

for j  = 1:6
    Data06reps(:,:,j) = mean(Data12reps(:,:, ...
                        index_to_start(j) : index_to_start(j) + 1),3);
end

end

