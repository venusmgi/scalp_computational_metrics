function epochIndices = findCleanEpochs(EEG_dur, fs, autoArts, epochLength, nToleratedArtSegment, durToleratedArtSegment)
%
% !!NOTE: Functionality for nToleratedArtSegment and durToleratedArtSegment
% has not been implemented yet!!!
%
% This function compares the duration of an EEG dataset to the times at
% which artifacts were detected, and it identifies clean epochs of data
% (no artifacts detected) of a desired duration. It has the option to
% "tolerate" some amount of artifact in the epoch, which can increase the
% amount of data available for analysis.
%
% INPUTS:
%   EEG_dur = duration of EEG dataset, in number of data points
%   fs = sampling frequency
%   autoArts = ouput of automatic artifact detection on the EEG dataset,
%       using the "get_automatedArtifacts_EEG" function
%   epochLength = desired duration of clean epochs, in seconds
%   nToleratedArtSegment = number of subsegments in the epoch that are
%       allowed to have artifacts
%   durToleratedArtSegment = duration of the epoch subsegment that is
%       allowed to contain artifacts
%
%   EXAMPLE: if epochLength=30, nToleratedArtSegment=2, and durToleratedArtSegment=1,
%   then each epoch of 30 seconds will be divided into 1-second
%   subsegments, and 2 of those subsegments can contain artifacts while
%   still deeming the epoch to be "clean." NOTE that this does NOT mean the
%   epoch can contain 2 artifacts, each 1-second long; it standardizes the
%   division of the epoch into 1-second subsegments first. This way, the
%   inputs can be selected to guarantee a number of clean subsegments
%   available for analysis of computational metrics.
%
% OUTPUT:
%   epochIndices = indices in the EEG data representing the start times of
%       clean epochs

% Create a binary vector the same length as the EEG; ones will represent artifacts
binArts = zeros(EEG_dur,1); 

% Create a vector of ones the same length as the epoch
binEpoch = ones(epochLength*fs,1);

% Loop through all artifacts, place ones at indices during artifacts
for i=1:size(autoArts.times,1)
    artInd = round(autoArts.times(i,:)*fs) + 1; % Note: artifacts are not necessarily in sequential order; shouldn't matter here
    binArts(artInd(1):artInd(2)) = 1;
end

%% Calculate indices of clean epochs

% % Example to test this section
% binArts = [1 0 0 0 0 0 0 0 0 1 1 1 0 0 0 1 1 1 1 1 1 0 0 0 0 0 1 1 0 0 0 1];
% binEpoch = [1 1 1 1];
% nToleratedArtSegment = 0;
% durToleratedArtSegment = 0;
% fs = 1;

% The result of this convolution indicates how many artifact points the
% epoch would overlap with, if the epoch started at that index
overlap = conv(binArts, binEpoch, 'valid');

maxOverlap = nToleratedArtSegment*durToleratedArtSegment*fs; % max # points overlapping with artifact
candidateStart = find(overlap <= maxOverlap); % potential start time for epoch

epochIndices = [];
ii = 1;
while ii<=length(candidateStart)
    epochIndices = [epochIndices candidateStart(ii)];  % save the index
    next_ii = find(candidateStart >= (candidateStart(ii)+length(binEpoch)), 1, 'first'); % the next eligible index can't overlap with the one we just saved
    if ~isempty(next_ii)
        ii = next_ii;  % if we find another eligible index, set ii to that value to it will get saved on the next iteration
    else
        ii = length(candidateStart)+1; % if there is no other eligible index, set ii to a high value so the loop is exited
    end
end

%% To add later
% If we want to add the variables nToleratedArtSegment and
% durToleratedArtSegment, we would next test every epoch identified above.
% Each epoch would be divided into segments of length
% durToleratedArtSegment, and we would test whether each sub-epoch was
% clean (sum = 0). If the number of clean sub-epochs is sufficient, then
% keep it, and if not, delete it from the epochIndices list.
% We would also need to save a vector or matrix indicating which sub-epochs
% are clean.