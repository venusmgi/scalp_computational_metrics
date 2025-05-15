function indices = Find_Clean_Indices(EEGDuration,fs,binArts,cleanEpochDuration)

% FIND_CLEAN_INDICES Identifies clean epochs in an EEG recording
%
% This function detects clean epochs in an EEG recording by identifying
% periods free of artifacts.
%
% Inputs:
%   EEGDuration - Duration of the EEG data in data points
%   fs - Sampling frequency in Hz
%   binArts - binary array indicating if an index is clean (0) or
%             artifactual (1)
%
%   cleanEpochDuration - Duration of desired clean epochs in seconds
%
% Output:
%   indices - Start indices of clean epochs



% Initialize variables for clean epoch detection
cleanEpochSamples = floor(cleanEpochDuration*fs);  % duration of clean epoch in samples
count = 0;  % counter used to measure duration of clean segment
indices = [];  % start indices of clean epochs

% Iterate through the binary artifact vector to find clean epochs
for i = 1:EEGDuration
    if binArts(i) == 0  % If it is a clean sample
        count = count + 1;  % add one to counter
        if count == cleanEpochSamples % if we have enough sequential clean samples
            indices(end+1,1) = i - cleanEpochSamples + 1;  % add starting index to vector
            count = 0;  % reset duration of clean segment to zero
        end
    else  % Artifact detected
        count = 0;  % reset duration of clean segment to zero
    end
end

return