function [nEpochs, startInd, stopInd] = Calc_Epoch_Indices(nSamp, fs, epochLength, startingIndices)
%
% This function calculates the number of epochs of a given length for an
% EEG record. It also checks to make sure none of the epochs extend beyond
% the length of the EEG.
% 
% INPUTS:
%   nSamp = number of samples (time) in EEG record
%   fs = sampling rate in Hz
%   epochLength = duration of epoch, in seconds
%   startingIndices = (optional) vector containing the starting indices for
%                     each epoch
%
% OUTPUTS:
%   nEpochs = number of epochs
%   startInd = vector containing starting index for each epoch
%   stopInd = vector containing ending index for each epoch


% Determine start and stop indices for each epoch
if nargin == 4
    assert(all(startingIndices> 0), 'Indices cannot be zero or negative; check the startingIndices values and make sure the smallest index is not smaller than 1.')
    nEpochs = length(startingIndices);  % number of epochs in data
    startInd = startingIndices; % starting index for each epoch
    stopInd = startingIndices+epochLength*fs - 1; % ending index for each epoch

else
    nEpochs = floor(nSamp./(epochLength*fs));  % number of epochs in data
    startInd = (0:nEpochs-1)*fs*epochLength + 1;  % starting index for each epoch
    stopInd = (1:nEpochs)*epochLength*fs;  % ending index for each epoch
end
if (~isempty(startInd) & (startInd(1) > fs*epochLength*10) ) %if the starting index is larger than 10 epochLength throw a warning
    warning('The starting clean index of this array is large, which may result in less clean epochs available.')
end

% Ensure indices are within bounds
oob = find(stopInd > nSamp);  % find epochs that are out of bounds (o.o.b)
nEpochs = nEpochs - length(oob);  % reduce number of epochs
startInd(oob) = [];  % remove start times for o.o.b. epochs
stopInd(oob) = [];   % remove stop times for o.o.b. epochs

