function [nEpochs, startInd, stopInd] = Calc_Epoch_Indices(nSamp, fs, epochLength, startingIndices)
%
% This function computes the starting and stopping sample indices for epochs of specified
% length within an EEG record. It also verifies that the calculated epochs fit within the
% total number of samples. If starting indices are provided, it bases the epoch division on
% these specified starting points.
% 
% Inputs:
%   nSamp - Total number of samples in the EEG record.
%   fs - Sampling rate of the EEG data in Hz.
%   epochLength - Desired length of each epoch in seconds.
%   startingIndices - Optional; a vector containing starting indices for each epoch. 
%                     If empty or not provided, the function will calculate indices based on epochLength.
%
% Outputs:
%   nEpochs - Total number of valid epochs detected or calculated.
%   startInd - Vector containing the start index for each valid epoch.
%   stopInd - Vector containing the stop index for each epoch.


% Check if startingIndices are provided, then use them to determine epoch indices
if nargin == 4
    % Error handling: Ensure that all starting indices are positive and within bounds
    assert(all(startingIndices> 0), 'All starting indices must be positive. Please check "startingIndices" and ensure no index is zero or negative.')
    nEpochs = length(startingIndices);  % number of epochs in data
    startInd = startingIndices; % starting index for each epoch
    stopInd = startingIndices+epochLength*fs - 1; % ending index for each epoch

else
    % Calculate equal-length epoch indices based on desired length
    nEpochs = floor(nSamp./(epochLength*fs));  % number of epochs in data
    startInd = (0:nEpochs-1)*fs*epochLength + 1;  % starting index for each epoch
    stopInd = (1:nEpochs)*epochLength*fs;  % ending index for each epoch
end

% Ensure indices are within bounds
oob = find(stopInd > nSamp);  % find epochs that are out of bounds (o.o.b)
if ~isempty(oob)
     warning(['Some epoch indices exceed the total sample length and are flagged as out of bounds.', ...
                 ' These indices have been removed.']);
end
nEpochs = nEpochs - length(oob);  % reduce number of epochs
startInd(oob) = [];  % remove start times for o.o.b. epochs
stopInd(oob) = [];   % remove stop times for o.o.b. epochs


