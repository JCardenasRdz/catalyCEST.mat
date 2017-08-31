function [ ImagingData, Bruker_Info] = load_data(message)
%
% LOAD_DATA : Wrapper to ask user to get the data needed.
% 
% SYNTAXIS
% ImagingData = LOAD_DATA(message)
% 
% INPUT
% message = string for the user to see

% display message 
uiwait(msgbox(message));


% get data after setting up pathway
experiment_directory=uigetdir(pwd); %Define folder pathname
[ImagingData, Bruker_Info] = Bruker_reader(experiment_directory);

ImagingData = squeeze(ImagingData);

end

