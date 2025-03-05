function [SEF,deltaDB,thetaDB,alphaDB,betaDB,broadDB] = Calc_SEF_SpectralPower_EEG(data,fs,epochLength,startingIndices)
% CALC_SEF_EEG_VENUS Calculates spectral edge frequency and band powers for EEG data
%
% Inputs:
%   data        - EEG data, channels x time
%   fs          - Sampling rate (Hz)
%   epochLength - Epoch length (seconds, default 5)

% Outputs:
%   SEF         - Spectral edge frequency value for each channel and epoch
%   deltaDB    - power of signgla in the delta band
%   bandPowers  - Structure containing power in different frequency bands
%
% Author: Venus Mostaghimi
% Date: 1.08.2025
% Version: 2.0


% Inputs:   data: EEG data, channels x time
%           fs: Sampling rate
%           epochLength: Epoch length (seconds, default 5)
% Outputs:  freq95: Spectral edge frequency value for each channel
%           powerMatFull: Matrix containing the Spectral edge frequency
%                         values for each channel up to 55 Hz
%
% Previously called: "calculateSEF_UCB.m"
% Code used in Smith et al. 2021
% Rachel J. Smith (Lopouratory 2019)
%
% Changelog
% V.1.0: Comment code, generalized code to fit generic EEG datsets
%        Remove 19 channel label setup from function
%        Combined calculateSEF_UCB and calculateSEF_UCB_full together
%        (both are identical but separately for some reason)
% V.1.1: Uploaded onto github (DH)
% V.1.2 : Venus added powerMat as the output, set fs to 200, removed
% resampling
% V.2.0 :Refactored for improved readability and performance (Venus Mostaghimi, 1.08.2025)



%Inpput Validation
validateattributes(fs, {'numeric'}, {'scalar','positive'}, 'Calc_SEF_EEG_Venus', 'fs', 2)
validateattributes(epochLength, {'numeric'}, {'scalar', 'positive'}, 'Calc_SEF_EEG_Venus','epochLength',3)

if nargin == 4
    nEpochs = length(startingIndices);
else
    % Initialize variables for outputs
    nSecs = floor(size(data,2)./fs);  % duration of data matrix, in seconds
    nEpochs = floor(nSecs./epochLength);  % number of epochs in data
end


nChans = size(data,1);  % number of EEG channels
SEF = nan(nEpochs,nChans);

% or define the variables one by one
[deltaDB,thetaDB,alphaDB,betaDB,broadDB] = deal(nan(nEpochs,nChans));


% Set 55 Hz cutoff for SEF calculation
cutoff = Find_Freq_Ind(55,epochLength);

% Set frequency ranges for EEG bands
deltaInd = Find_Freq_Ind([1 4],epochLength);
thetaInd = Find_Freq_Ind([4 8],epochLength);
alphaInd = Find_Freq_Ind([8 13],epochLength);
betaInd = Find_Freq_Ind([13 30],epochLength);


%% I can also define this as a struct:
bandPowers = struct('delta', deltaDB, 'theta', thetaDB, 'alpha', alphaDB, 'beta', betaDB, 'broad', broadDB);

% For each epoch, calculate the power spectrum
for epochInd = 1:nEpochs
    %% Select epoch
    if nargin == 4
        startInd = startingIndices(epochInd);
        stopInd = startingIndices(epochInd)+epochLength*fs-1;
    else
        startInd = (epochInd-1)*(epochLength*fs)+1;
        stopInd = (epochInd*epochLength*fs);
    end

    % Ensure indices are within bounds
    if stopInd > size(data, 2)
        % stopInd = size(data, 2);
        break; % Exit the loop, because I don't want info from an epoch size that is shorter than others
    end


    eegEpoch = data(:,startInd:stopInd)'; % transpose data so each channel is in a column

    % Apply windowing function to epoch of EEG data
    %% VENUS CHECK THIS, I HAVEN'T DEBUGGED IT
    win = hann(stopInd-startInd+1);  % create hanning window based on number of time points in eegEpoch
    eegEpochWin = repmat(win,1,nChans).*eegEpoch;  % apply window to EEG in each channel

    % Calculate FFT and EEG power
    fftVals = fft(eegEpochWin);

    % Calculate power; dimension is # frequencies x # channels
    powerVals = abs(fftVals).^2; 
    powerVals = powerVals(1:length(powerVals)/2+1, :);  % remove negative frequencies (two-sided power to one-sided power)
    powerVals(2:end-1) = 2*powerVals(2:end-1); % doubling the power except for DC and Nyquist

    % Calculate SEF
    SEF(epochInd,:) = Calc_SEF(powerVals,fs,cutoff);  % dimension is 1 x channels

    % Calculate power in each frequency band

    deltaDB(epochInd,:) = Calc_Total_Power(powerVals, deltaInd, fs, win);
    thetaDB(epochInd,:) = Calc_Total_Power(powerVals, thetaInd, fs, win);
    alphaDB(epochInd,:) = Calc_Total_Power(powerVals, alphaInd, fs, win);
    betaDB(epochInd,:)  = Calc_Total_Power(powerVals, betaInd, fs, win);
    broadDB(epochInd,:) = Calc_Total_Power(powerVals, 'all', fs, win);


end

return

%% do we want to turn SEF ro decible as well?
% Function to calculate SEF for one epoch
function SEF = Calc_SEF(powerVals,fs, cutoff)
% CALC_SEF Calculates the Spectral Edge Frequency (SEF) for one epoch
%
% Inputs:
%   powerVals - Matrix of FFT values, frequencies x channels
%   cutoff  - Index corresponding to the maximum frequency to consider (e.g., 55 Hz)
%   fs  - Sampling frequency (e.g. 200Hz)
%
% Output:
%   SEF     - Spectral Edge Frequency for each channel (95th percentile)

powerSpec = powerVals(1:cutoff,:);  % frequencies (up to 55 Hz) x channels

% Define frequency range from 0 to 55 Hz
frequencies = linspace(0, fs/2, length(powerVals)); %% DR.LOPOUR, shouldn't i have up to 55 Hz instead of fs/2?

% Calculate the cumulative sum of the power spectrum
cumulativePower = cumsum(powerSpec);

% Normalize the cumulative power to create a cumulative distribution function (CDF)
cumulativePowerCDF = cumulativePower / cumulativePower(end);

% Find the index where the CDF reaches or exceeds the 95th percentile
percentileIndex = find(cumulativePowerCDF >= 0.95, 1);

% Find the frequency corresponding to this index
SEF = frequencies(percentileIndex);


return

%% I can also define it as below:
% Calc_SEF = @(powerVals, cutoff) prctile(powerVals(1:cutoff,:), 95, 1);


% Function to find the index for a specific frequency value associated with the fft
% result
function fInd = Find_Freq_Ind(freq, epochLength)
% FIND_FREQ_IND Finds the index for a specific frequency value in FFT results
%
% Inputs:
%   freq        - Desired frequency (Hz)
%   epochLength - Length of the epoch (seconds)
%
% Output:
%   fInd        - Index corresponding to the input frequency in FFT results
%
fInd = floor(freq * epochLength);
% Alternative1 calculation
% fInd = round(freq * N / fs) + 1;
% Alternative2 calculation for readability:
% nyqF = fs/2;
% fInd = floor((freq * nyqF * epochLength) / nyqF)
return

%% I can also define it as below:
% Find_Freq_Ind = @(freq, epochLength) floor(freq*epochLength)

function powerDecibel = Calc_Total_Power(powerVals, freqInd, fs, win)
% CALC_TOTAL_POWER Calculates the total power in a specified frequency band
%
% Inputs:
%   powerVals     - Matrix of FFT values, frequencies x channels
%   freqInd     - Frequency indices to consider. Can be:
%                 - 'all' for all frequencies
%                 - [start_index, end_index] for a specific frequency range
%   epochLength - Length of the epoch (seconds)
%   fs          - Sampling frequency (Hz)
%   win         - Window function used in the FFT calculation
%
% Output:
%   powerDecibel - Total power in the specified frequency band (dB)
%
% Note: This function calculates the power spectral density and then
% integrates it over the specified frequency range to get total power.


% Extract power values for the desired frequency band
if strcmp(freqInd, 'all')
    powerMat = powerVals(:, :);
else
    powerMat = powerVals(freqInd(1):freqInd(2), :);  % dimension is frequencies x channels
end

% Convert power to power spectral density
S = sum(win.^2);  % scaling factor based on hanning window function
psdEEG = powerMat /( fs * S);  % Divide by fs and multiply by scaling factor to obtain power spectral density (power/Hz)

% Calculate total power in the desired frequency band and convert to
% decibels
powerLinear = trapz(psdEEG,1);  % calculate the area; 
powerDecibel = 10*log10(powerLinear);  % convert total power to decibels
return
