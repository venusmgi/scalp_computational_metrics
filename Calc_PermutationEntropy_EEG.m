function [permEnt] = Calc_PermutationEntropy_EEG(data,fs,order,delay,epochLength,startingIndices)


% Calculates the Permutation Entropy for all EEG channels and epochs.
%
% Inputs:
%   data       : Filtered clean EEG data (channels x time points).
%   fs             : Sampling frequency of the EEG data.
%   order          : Order of permutation entropy (default 5).
%   delay          : Delay time of permutation entropy (default 1).
%   epochLength    : Length of each epoch in seconds.
%   startingIndices: (Optional) Starting indices for epochs.
%
% Outputs:
%   permEnt        : Permutation entropy for each channel and epoch.

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

nChanns = size(data,1);
if nargin == 6
    nEpochs = length(startingIndices);
else
    nEpochs = floor(size(data,2)./(epochLength*fs));  % number of epochs in data
end

% Initialize permutation entropy matrix
permEnt = nan(nEpochs, nChanns);

    for chanId = 1: nChanns
        for epochId = 1:nEpochs
            if nargin == 6
                startInd = startingIndices(epochId);
                stopInd = startInd + epochLength * fs - 1;
            else
                startInd = (epochId - 1) * epochLength * fs + 1;
                stopInd = epochId * epochLength * fs;
            end

                        % Ensure indices are within bounds
            if stopInd > size(data, 2)
                % stopInd = size(data, 2);
                break; % Exit the loop, because I don't want info from an epoch size that is shorter than others
            end

            [permEnt(epochId,chanId),~] = Pec(data(chanId,startInd:stopInd),order,delay);
        end
    end
    


return
