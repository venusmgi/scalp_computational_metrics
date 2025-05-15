function indices = Find_Clean_Indices(EEGDuration,fs,autoArts,cleanEpochDuration)

% FIND_CLEAN_INDICES Identifies clean epochs in an EEG recording
%
% This function detects clean epochs in an EEG recording by identifying
% periods free of artifacts.
%
% Inputs:
%   EEGDuration - Duration of the EEG data in data points
%   fs - Sampling frequency in Hz
%   autoArts - artifact times; this is the output of the
%              get_automatedArtifacts_EEG function
%   cleanEpochDuration - Duration of desired clean epochs in seconds
%
% Output:
%   indices - Start indices of clean epochs

% Create a binary vector representing artifacts (1 = artifact, 0 = clean)
binArts = zeros(EEGDuration, 1);
for i = 1:size(autoArts.times,1)
    artStart = round(autoArts.times(i, 1)*fs)+1;
    artEnd = min(EEGDuration, round(autoArts.times(i, 2) * fs));
    binArts(artStart:artEnd) = 1; % Note: artifacts are not necessarily in sequential order; shouldn't matter here
end

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