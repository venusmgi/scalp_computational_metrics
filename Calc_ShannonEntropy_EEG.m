function entropy = Calc_ShannonEntropy_EEG(data,fs,aveOptimal,epochLength,startingIndices)

% Calculates the Shannon entropy value for each channel.
% 
% Inputs:
%   data       : Filtered clean EEG data (channels x time points).
%   aveOptimal     : Optimal bin count to use (default 390). Use the 
%                    Freedman-Diaconis rule on the signal to calculate:
%                    nBins = (max(x) - min(x)) / (2 * IQR(x) * length(x)^(-1/3)).
%   epochLength    : Length of each epoch in seconds.
%   startingIndices: (Optional) Starting indices for epochs.
%
% Output:
%   entropy        : Shannon entropy value for each channel.
%
% Previously called: "calculateEntropy.m"
% Code used in Smith et al. 2021
% Rachel J. Smith (Lopouratory 2019)
%
% V.1.0: Derek Hu Comment code, generalized code to fit generic EEG datasets.
% V.1.1: Derek Hu Uploaded onto GitHub.
% V.2.0: Venus added option to apply entropy calculation for epochs.


% Number of channels
    nChanns = size(data, 1);

    % Determine number of epochs
    if nargin == 5
        nEpochs = length(startingIndices);
    else
        nEpochs = floor(size(data,2)./(epochLength*fs));  % number of epochs in data
    end


    % Initialize entropy matrix
    entropy = nan(nEpochs, nChanns);

     % Calculate entropy for each channel and epoch
    for chan = 1:nChanns
        for epochInd = 1:nEpochs
            if nargin == 5
                startInd = startingIndices(epochInd);
                stopInd = startInd + epochLength * fs - 1;
            else
                startInd = (epochInd - 1) * epochLength * fs + 1;
                stopInd = epochInd * epochLength * fs;
            end
            
            % Ensure indices are within bounds
            if stopInd > size(data, 2)
                % stopInd = size(data, 2);
                break; % Exit the loop, because I don't want info from an epoch size that is shorter than others
            end
            
            % Calculate histogram and probabilities
            [N, ~] = histcounts(round(data(chan, startInd:stopInd),4), aveOptimal);
            probs = N ./ sum(N);
            entVal = -sum(probs .* log2(probs + eps));
            
            % Store entropy value
            entropy(epochInd, chan) = entVal;
        end
    end
return