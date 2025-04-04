function ampMatrix = Calc_Range_EEG(data, fs, epochLength, startingIndices)
% Outputs a matrix of all amplitude values for each epoch and for all
% channels in the input data, calculated using the range (peak-to-peak).
% If artifactual epochs need to be removed, they can be identified using
% get_automatedArtifacts_EEG and the artifactual amplitude values can be
% removed later.
%
% Inputs:
%   data        - EEG data, channels x time points
%   fs          - Sampling rate
%   epochLength - Epoch length (seconds; suggested value=1)
%   startingIndices - (Optional) Vector of starting indices for each epoch. If not provided,
%                     epochs are determined based on the epoch length.

% Outputs:
%    ampMatrix: matrix of all amplitude values for each one second window
%               (nEpochs x nChannels)
% Notes:
% - The function assumes the input EEG data is already preprocessed
% - The function will automatically handle the calculation of epochs based on
%   the provided epoch length or starting indices. If starting indices are provided,
%   they must be consistent with the epoch length to avoid overlap or gaps.
%   If the epoch length is larger/smaller than the distance between each
%   index in starting indices, then you will have over lap/gaps.

%
% Previously called: "calculate_Amplitude_full.m"
% Code used in Smith et al. 2021
% Rachel J. Smith (Lopouratory 2019)
%
% V.1.0: Comment code, generalized code to fit generic EEG datsets (DH)
%        Remove 19 channel label setup from function
% V.1.1: Uploaded onto github (DH)
% V2: Updated based on lab code review on 12/18/24

% IDEAS FOR SUBSEQUENT REVISION:
% - transpose data matrix in main code for faster processing?

nChan = size(data,1);  % number of EEG channels
nSamp = size(data,2);  % number of samples (time)

if nargin == 4
    [nEpochs, startInd, stopInd] = Calc_Epoch_Indices(nSamp, fs, epochLength, startingIndices);
else
    [nEpochs, startInd, stopInd] = Calc_Epoch_Indices(nSamp, fs, epochLength);
end

% Check that the data matrix is in the correct orientation
assert(nChan < nSamp, 'Warning: Number of channels is greater than number of time samples! Data matrix may need to be transposed.')

% Check that the data contains at least 1 second of data
assert(nSamp >= fs, 'Warning: Data matrix contains less than one second of data.')

% Initialize matrix for amplitude values, epochs x channels
ampMatrix = nan(nEpochs,nChan);

% Loop through all epochs and calculate the range for each one
for win=1:nEpochs
    eegWin = data(:,startInd(win):stopInd(win));  % select epoch of EEG
    ampMatrix(win,:) = (max(eegWin,[],2) - min(eegWin,[],2))'; % Calculate range for this window of data
end


