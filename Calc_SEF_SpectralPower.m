function [SEF,deltaDB,thetaDB,alphaDB,betaDB,broadDB] = Calc_SEF_SpectralPower_EEG(data,fs,epochLength,startingIndices)
% CALC_SEF_SPECTRALPOWER_EEG Calculates spectral edge frequency and band powers for EEG data
%
% Inputs:
%   data        - EEG data, channels x time
%   fs          - Sampling rate (Hz)
%   epochLength - Epoch length (seconds, suggested value=5)
%   startingIndices - (Optional) Vector of starting indices for each epoch. If not provided,
%                     epochs are determined based on the epoch length.
%
% Outputs:
%     *** NOTE: All outputs have dimensions nEpochs x nChan
%   SEF           - Spectral Edge Frequency for each channel and epoch, representing the frequency
%                   below which 95% of the power is contained.
%   deltaDB       - Power of the signal in the delta band (1-4 Hz) expressed in decibels.
%   thetaDB       - Power of the signal in the theta band (4-8 Hz) expressed in decibels.
%   alphaDB       - Power of the signal in the alpha band (8-13 Hz) expressed in decibels.
%   betaDB        - Power of the signal in the beta band (13-30 Hz) expressed in decibels.
%   broadDB       - Power of the signal in the broadband (full spectrum) expressed in decibels.
%
% Notes:
% - The function assumes the input EEG data is already preprocessed and ready
%   for spectral analysis.
% - The function will automatically handle the calculation of epochs based on
%   the provided epoch length or starting indices. If starting indices are provided,
%   they must be consistent with the epoch length to avoid overlap or gaps.
%   If the epoch length is larger/smaller than the distance between each
%   index in starting indices, then you will have over lap/gaps.
% - The SEF calculation uses a cutoff frequency of 55 Hz, which can be adjusted if necessary.

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

%Input Validation
validateattributes(fs, {'numeric'}, {'scalar','positive'}, 'Calc_SEF_SpectralPower_EEG', 'fs', 2)
validateattributes(epochLength, {'numeric'}, {'scalar', 'positive'}, 'Calc_SEF_EEG_Venus','epochLength',3)

% Parameters
SEFCutoff = 55; % frequency cutoff for SEF, in Hz

% Number of EEG channels and samples in input data
nChan = size(data,1);
nSamp = size(data,2);

% Get epoch start and stop indices
if nargin == 4
    [nEpochs, startInd, stopInd] = Calc_Epoch_Indices(nSamp, fs, epochLength, startingIndices);
else
    [nEpochs, startInd, stopInd] = Calc_Epoch_Indices(nSamp, fs, epochLength);
end

% Initialize variables for outputs
SEF = nan(nEpochs,nChan);
[deltaDB,thetaDB,alphaDB,betaDB,broadDB] = deal(nan(nEpochs,nChan));

% Set 55 Hz cutoff for SEF calculation
cutoffInd = Find_Freq_Ind(SEFCutoff,epochLength);

% Set frequency ranges for EEG bands
deltaInd = Find_Freq_Ind([1 4],epochLength);
thetaInd = Find_Freq_Ind([4 8],epochLength);
alphaInd = Find_Freq_Ind([8 13],epochLength);
betaInd = Find_Freq_Ind([13 30],epochLength);

% Create hanning window based on number of time points in EEG epoch; column vec
win = hann(epochLength*fs);  

% For each epoch and each channel, calculate the power spectrum
for epochInd = 1:nEpochs
    % Extract the EEG epoch for all channels
    eegEpoch = data(:,startInd(epochInd):stopInd(epochInd))'; % transpose data so each channel is in a column
    % Apply windowing function to epoch of EEG data
    eegEpochWin = repmat(win,1,nChan).*eegEpoch;  % apply window to EEG in each channel

    % Calculate FFT and EEG power
    fftVals = fft(eegEpochWin);

    % Calculate power; dimension is # frequencies x # channels
    powerVals = abs(fftVals).^2;
    powerVals = powerVals(1:length(powerVals)/2+1, :);  % remove negative frequencies (two-sided power to one-sided power)
    powerVals(2:end-1,:) = 2*powerVals(2:end-1,:); % doubling the power except for DC and Nyquist

    % Calculate SEF
    SEF(epochInd,:) = Calc_SEF(powerVals, fs, cutoffInd);

    % Calculate power in each frequency band
    deltaDB(epochInd,:) = Calc_Band_Power(powerVals, deltaInd, fs, win);
    thetaDB(epochInd,:) = Calc_Band_Power(powerVals, thetaInd, fs, win);
    alphaDB(epochInd,:) = Calc_Band_Power(powerVals, alphaInd, fs, win);
    betaDB(epochInd,:)  = Calc_Band_Power(powerVals, betaInd, fs, win);
    broadDB(epochInd,:) = Calc_Band_Power(powerVals, 'all', fs, win);
end

return

%% Function to calculate SEF for one epoch
function SEF = Calc_SEF(powerVals, fs, cutoffInd)
% CALC_SEF Calculates the Spectral Edge Frequency (SEF) for one epoch
%
% Inputs:
%   powerVals - Matrix of FFT values for all positive frequencies; frequencies x channels
%   cutoffInd  - Index corresponding to the maximum frequency to consider (e.g., 55 Hz)
%   fs  - Sampling frequency, in Hz (e.g. 200Hz)
%
% Output:
%   SEF  - Spectral Edge Frequency for each channel (95th percentile)

% Extract all frequencies up to cutoff freq; dimension is frequencies x channels
powerSpec = powerVals(1:cutoffInd,:);

% Define frequency range for powerVals matrix (negative frequencies have
% already been removed)
frequencies = linspace(0, fs/2, size(powerVals,1));

% Calculate the cumulative sum of the power spectrum in each column
cumulativePower = cumsum(powerSpec);

% Normalize the cumulative power to create a cumulative distribution function (CDF)
cumulativePowerCDF = cumulativePower./repmat(cumulativePower(end,:),size(powerSpec,1),1);

% Find the index where the CDF reaches or exceeds the 95th percentile
[~,percentileIndex] = min(abs(cumulativePowerCDF-0.95));  % output is row vector

% Find the frequency corresponding to the index for each channel
SEF = frequencies(percentileIndex);  % row vector, 1 x nChan

return

%% Function to find index for a specific frequency value associated with fft result
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

% Alternative1 calculation
fInd = round(freq * epochLength)+1; %This method has the smallest error
% Alternative2 calculation for readability:
% nyqF = fs/2;
% fInd = floor((freq * nyqF * epochLength) / nyqF)
% Alternative 3
% Find_Freq_Ind = @(freq, epochLength) floor(freq*epochLength)
return

%% Function to calculate total power in specified frequency band
function powerDecibel = Calc_Band_Power(powerVals, freqInd, fs, win)
% CALC_BAND_SPECIFIC_TOTAL_POWER Calculates the total power in a specified frequency band
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

% Frequency vector for all PowerVals frequencies
allFreq = linspace(0,fs/2,size(powerVals,1));

% Extract power values for the desired frequency band
if strcmp(freqInd, 'all')
    powerMat = powerVals;
    bandFreq = allFreq;
else
    powerMat = powerVals(freqInd(1):freqInd(2), :);  % dimension is frequencies x channels
    bandFreq = allFreq(freqInd(1):freqInd(2));  % frequencies for specified band
end

% Convert power to power spectral density
S = sum(win.^2);  % scaling factor based on hanning window function (scalar)
psdEEG = powerMat/(fs * S);  % Divide by fs and scaling factor to obtain power spectral density (power/Hz)

% Calculate total power in the desired frequency band and convert to decibels
powerLinear = trapz(bandFreq, psdEEG,1);  % calculate the area
powerDecibel = 10*log10(powerLinear);  % convert total power to decibels; dimension is 1 x nChan

return
