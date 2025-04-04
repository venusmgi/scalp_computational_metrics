function entropy = Calc_ShannonEntropy_EEG(data,fs,nBins,epochLength,startingIndices)

% Calculates the Shannon entropy value for each channel.
%
% Inputs:
%   data       : Filtered clean EEG data (channels x time points).
%   nBins      : Optimal bin count to use (default 390). Use the
%                    Freedman-Diaconis rule on the signal to calculate:
%                    nBins = (max(x) - min(x)) / (2 * IQR(x) * length(x)^(-1/3)).
%   epochLength    : Length of each epoch in seconds.
%   startingIndices: (Optional) Starting indices for epochs.
%
% Output:
%   entropy  : Shannon entropy values; dimension is nEpochs x nChan
%
% Previously called: "calculateEntropy.m"
% Code used in Smith et al. 2021
% Rachel J. Smith (Lopouratory 2019)
%
% V.1.0: Derek Hu Comment code, generalized code to fit generic EEG datasets.
% V.1.1: Derek Hu Uploaded onto GitHub.
% V.2.0: Venus added option to apply entropy calculation for epochs.

% Number of channels and samples in EEG record
nChan = size(data, 1);
nSamp = size(data,2);

% Determine start and stop indices for each epoch
if nargin == 5
    [nEpochs, startInd, stopInd] = Calc_Epoch_Indices(nSamp, fs, epochLength, startingIndices);
else
    [nEpochs, startInd, stopInd] = Calc_Epoch_Indices(nSamp, fs, epochLength);
end

% Initialize entropy matrix
entropy = nan(nEpochs, nChan);

% Calculate entropy for each channel and epoch
for chanId = 1:nChan
    for epochInd = 1:nEpochs
        % Extract EEG epoch
        epochEEG = data(chanId, startInd(epochInd):stopInd(epochInd));

        % Calculate histogram and probabilities
        [N, ~] = histcounts(epochEEG, nBins);
        probs = N ./ sum(N);
        epochEntropy = -sum(probs .* log2(probs + eps)); % entropy for one epoch and one channel

        % Store entropy value
        entropy(epochInd, chanId) = epochEntropy;
    end
end

return