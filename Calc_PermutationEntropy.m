function permEnt = Calc_PermutationEntropy(data,fs,order,delay,epochLength,startingIndices)
% Calculates the Permutation Entropy for all EEG channels and epochs.
%
% Inputs:
%   data     : Filtered clean EEG data (channels x time points).
%   fs       : Sampling frequency of the EEG data.
%   order    : Order of permutation entropy (default 5).
%   delay    : Delay time of permutation entropy (default 1).
%   epochLength    : Length of each epoch in seconds.
%   startingIndices: (Optional) Starting indices for epochs.
%
% Outputs:
%   permEnt  : Permutation entropy for each channel and epoch;
%               dimension is nEpoch x nChan

%
% Previously called: "calculate_PermEntropy_Derek.m"
% Revised from code used in Smith et al. 2021
% Uses the pec function (Ouyang et al.)
% Derek Hu (Lopouratory, 2021)
%
% V.1.0: Comment code, generalized code to fit generic EEG datasets.
% V.1.1: Uploaded onto GitHub.
% V.2.0 : Venus added option to apply entropy calculation for epochs.
% Removed new_size from the output

% Number of channels and samples in EEG record
nChan = size(data,1);
nSamp = size(data,2);

% Determine start and stop indices for each epoch
if nargin == 6
    [nEpochs, startInd, stopInd] = Calc_Epoch_Indices(nSamp, fs, epochLength, startingIndices);
else
    [nEpochs, startInd, stopInd] = Calc_Epoch_Indices(nSamp, fs, epochLength);
end

% Initialize permutation entropy matrix
permEnt = nan(nEpochs, nChan);

% Loop through each channel and each epoch
for chanId = 1:nChan
    for epochId = 1:nEpochs
        % Extract EEG epoch
        epochEEG = data(chanId,startInd(epochId):stopInd(epochId));

        % Calculate permutation entropy
        [permEnt(epochId,chanId),~] = Pec(epochEEG,order,delay);
    end
end

return
