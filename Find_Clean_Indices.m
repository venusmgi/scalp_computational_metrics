function [indices] = Find_Clean_Indices(EEG_record,fs, clean_epoch_duration_in_second)

% FIND_CLEAN_INDICES Identifies clean epochs in an EEG recording
%
% This function detects clean epochs in an EEG recording by identifying
% periods free of artifacts.
%
% Inputs:
%   EEG_record - The EEG data (channels x time points)
%   fs - Sampling frequency in Hz
%   clean_epoch_duration_in_second - Duration of desired clean epochs in seconds
%
% Output:
%   indices - Start indices of clean epochs

% Parameters for artifact detection
stdAbove = 7.5;        % Standard deviation threshold for artifact detection
buffer = 0.9;          % Buffer around detected artifacts (in seconds)
channelsInvolved = 1;  % Number of channels that must exceed threshold for artifact detection

% Get the total duration of the EEG recording
EEG_dur = size(EEG_record, 2);

% Detect artifacts using the automated artifact detection function
autoArts = get_automatedArtifacts_EEG(EEG_record, fs, stdAbove, buffer, channelsInvolved);

% Create a binary vector representing artifacts (1 = artifact, 0 = clean)
binArts = zeros(EEG_dur, 1);
artInd = nan(size(autoArts.times,1),2);
for i = 1:size(autoArts.times,1)
    artStart = max(1, round(autoArts.times(i, 1) * fs));
    artEnd = min(EEG_dur, round(autoArts.times(i, 2) * fs));
    artInd(i,:) = [artStart,artEnd]; % Note: artifacts are not necessarily in sequential order; shouldn't matter here
    binArts(artStart:artEnd) = 1;
end

% Initialize variables for clean epoch detection
clean_epoch_samples = floor(clean_epoch_duration_in_second * fs);
count = 0;
indices = [];

% Iterate through the binary artifact vector to find clean epochs
for i = 1:EEG_dur
    if binArts(i) == 0  % Clean sample
        count = count + 1;
        if count == clean_epoch_samples %f mod(count, clean_epoch_samples)==0
            indices(end+1) = i - clean_epoch_samples + 1;
            count = 0;
        end
    else  % Artifact detected
        count = 0;
    end
end

end