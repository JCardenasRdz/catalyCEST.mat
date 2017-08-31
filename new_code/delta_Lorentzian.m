function [ delta_Lorentzian ] = delta_Lorentzian( parameters, PPM_matrix)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes hereX
%  delta_Lorentzian = Lsum_post - Lsum_pre;

npools = 2;

% pre
ppm_pre =  PPM_matrix(:,1);
parameters_pre =  parameters(1:npools:end);
Lsum_pre =lorentzian(parameters_pre,ppm_pre);

% post
ppm_post = PPM_matrix(:,2);
parameters_post = parameters(2:npools:end);
Lsum_post =lorentzian(parameters_post,ppm_post);

delta_Lorentzian = Lsum_post - Lsum_pre;

end

