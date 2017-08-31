%function [] = catalyCEST_MRI_invivoMAPS(  )
% ACIDOCEST_MRI_INVIVO performs acidoCEST MRI analysis assuming the
% following:
% 1) Isovue was used as the contrast agent at 3.0uT
% 2)
% 3)

%% Preinjection
% load  data
message = 'Load pre-injection data';
[ imgs_pre_scan, Bruker_Info] = load_data(message);
offsets_pre = Bruker_Info.cest_array   ; 

%% reshape and average     
[ imgs_pre_scan ] = reshape_and_average( imgs_pre_scan, offsets_pre );
%% Post-injection        
% load post-injection data
message = 'Load post-injection data';
[ imgs_post_scan, Bruker_Info] = load_data(message);
offsets_post = Bruker_Info.cest_array   ; 
     
% reshape and average     
[ imgs_post_scan_avg ] = reshape_and_average( imgs_post_scan, offsets_post );

%% select ROIs
[ mask_tumor, mask_noise] = select_ROIS( imgs_post_scan);

%% Control 
figure(1);
plot(offsets_pre, avgroi(imgs_pre_scan,mask_tumor)); hold all
plot(offsets_post + 0, avgroi(imgs_post_scan_avg,mask_tumor))
legend({'Pre-injection','Post'})

uiwait(msgbox('Does the average data looks good to you?'));
close all
%% reshape postscan
[rows,cols,num_images] = size(imgs_post_scan);
num_sat_freqs = length(offsets_post);
repetitions = num_images / num_sat_freqs;

imgs_post_scan = reshape(imgs_post_scan,[rows,cols,num_sat_freqs,repetitions]);

%% Run Lorentzian fitting
% Parameters for fitting
X0(:,1)=[0,0,2,2,4,9.3]'; %Lower bound guesses for parameters (amplitudes,widths,offset)
X0(:,2)=[0.02,0.02,2.5,2.5,4.8,9.8]'; %Initial guesses for parameters
X0(:,3)=[.15,.15,3,3,5.5,13.0]'; %Upper bound guesses for parameters

Method.Npools=2;      % Number of pools
Method.range=[3, num_sat_freqs-2]; %First and last point of the Zspectra to be analyzed
Method.x0=X0;         % Allocate matrix into stru


% allocate output
CEST1_cube= zeros(rows,cols,repetitions);
CEST2_cube = zeros(rows,cols,repetitions);
RSQ_cube =  zeros(rows,cols,repetitions);

ppm = offsets_post/ 300;

tic
for j= 1:repetitions
    
% pass each repetition
[~,MAPS,indices]=cestMRI_voxel(imgs_pre_scan,imgs_post_scan(:,:,:,j), ...
                                                ppm,mask_tumor,Method);

RSQ_cube(:,:,j) = MAPS.rsq;
CEST1_cube(:,:,j) = MAPS.cest(:,:,1);
CEST2_cube(:,:,j) = MAPS.cest(:,:,2);
end
toc


%% Reaction Coordinate
cest_9ppm_corrected = CEST2_cube;
cest_9ppm_corrected(CEST2_cube < 1E-3)  = NaN;

% create colormap to get black at first number in range
activity_colormap = jet(64);
activity_colormap(1,:) = 0;

% activity
RxCoordinate = 1 - CEST1_cube ./ cest_9ppm_corrected;

%% Load Rate
message = 'Load anatomical referencedata';
[ RARE, ~] = load_data(message);

%%

