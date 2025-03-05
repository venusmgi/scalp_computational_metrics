function [filteredEEG,startingIndCleanEpoch] = Rereference_Filter_DetectArts(EEG_record,EEG_hdr,rereferenceMethod,filterType,cleanEpochLength,artifactDetection)
% REREFERENCE_FILTER_DETECTARTS Re-references, filters, and detects clean epochs in EEG data
%
% This function processes EEG data by re-referencing, filtering, and
% detecting clean epochs across all channels.
%
% Inputs:
%   EEG_record - The EEG data (channels x time points)
%   EEG_hdr - Header information containing sampling frequency
%   rereferenceMethod - Method for re-referencing ('EAR' or 'CAR')
%   filterType - Type of filter to apply
%   cleanEpochLength - Length of desired clean epochs in seconds
%   artifactDetection - Boolean flag to perform artifact detection
%
% Outputs:
%   filteredEEG - The filtered EEG data
%   startingIndCleanEpoch - Start indices of clean epochs
%   cleanEEG - Clean EEG segments (channels x time points)

% Re-reference the EEG data
if strcmp(rereferenceMethod, 'EAR')
    %finding the location of the ear or mastiod channeld
    channel = 1:length(EEG_hdr.label);
    A1 = channel(strcmp(EEG_hdr.label ,'A1') | strcmp(EEG_hdr.label ,'M1'));
    A2 = channel(strcmp(EEG_hdr.label ,'A2') | strcmp(EEG_hdr.label ,'M2'));
    % Finding linked ear values
    linkedEar = (EEG_record(A1,:) + EEG_record(A2,:)) / 2; 
    EEG_reref = EEG_record - repmat(linkedEar, size(EEG_record, 1), 1); % Re-referencing the signal
elseif strcmp(rereferenceMethod, 'CAR')
    CAR = mean(EEG_record(1:19, :), 1); % Average over main 19 electrodes
    EEG_reref = EEG_record - repmat(CAR, size(EEG_record, 1), 1); % Re-referencing the signal
else
    error('Invalid re-reference method. Choose either "EAR" or "CAR".');
end

% Get the sampling frequency
fs = EEG_hdr.frequency(1);

% Filter the re-referenced signal
filter = Pick_Filter(filterType, fs);
filteredEEG = filtfilt(filter, 1, EEG_reref')';

if artifactDetection
 %finding clean epochs
 startingIndCleanEpoch = Find_Clean_Indices(filteredEEG,fs, cleanEpochLength)';
else
    startingIndCleanEpoch = [];
end
 % cleanEpochs = findCleanEpochs(filteredEEG, fs, epochLength, nToleratedArtSegment, durToleratedArtSegment);

 %% Extract clean EEG segments across all channels (this paert takes a long time, find a better way)
% 
%  % Calculate the number of samples per clean epoch
%  epochSamples = round(cleanEpochLength * fs);
% 
%  % Estimate the total number of clean samples
%  totalCleanSamples = epochSamples * length(startingIndCleanEpoch);
% 
%  % Preallocate the cleanEEG matrix
%  cleanEEG = nan(size(EEG_record,1), totalCleanSamples);
% 
%  currentPos = 1
%  for i = 1:length(startingIndCleanEpoch)
%      startInd = min((startingIndCleanEpoch(i)-1)*fs +1,1);
%      endInd = min(startingIndCleanEpoch(i)*fs, length(filteredEEG));
%        % Ensure indices are within bounds
%     if endInd <= size(filteredEEG, 2)
%          % Calculate the number of samples in the current epoch
%         numSamples = endInd - startInd + 1;
% 
%         cleanEEG(:, currentPos: currentPos+numSamples-1 ) =  filteredEEG(:, startInd:endInd);
%         currentPos = currentPos+numSamples;
%     end
% 
%  end
% 
%  % Trim the cleanEEG matrix to remove any unused preallocated space
% cleanEEG = cleanEEG(:, 1:currentPos - 1);

return