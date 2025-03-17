function [SEF,deltaDB,thetaDB,alphaDB,betaDB,broadDB] = Calc_SEF_SpectralPower_EEG(data,fs,epochLength,startingIndices)
% CALC_SEF_SPECTRALPOWER_EEG Calculates spectral edge frequency and band powers for EEG data
%
% Inputs:
%   data        - EEG data, channels x time
%   fs          - Sampling rate (Hz)
%   epochLength - Epoch length (seconds, default 5)
%   startingIndices - (Optional) Vector of starting indices for each epoch. If not provided,
%                     epochs are determined based on the epoch length.
%
% Outputs:
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

%vectorizing startInd and stopInd
if nargin == 4
    nEpochs = length(startingIndices);
    startInd = startingIndices;
    stopInd = startingIndices+fs*epochLength-1;
    %% QUESTION
    % should I remove this out of the if statement and regardless
    % of whatever it starting indecies if just put windowLength = epochLength*fs;?

    windowLength = unique(stopInd-startInd+1);
else

    nSecs = floor(size(data,2)./fs);  % duration of data matrix, in seconds
    nEpochs = floor(nSecs./epochLength);  % number of epochs in data
    startInd = (0:nEpochs-1)*epochLength*fs + 1;
    stopInd = (1:nEpochs)*epochLength*fs;
    windowLength = epochLength*fs;
end

% Ensure indices are within bounds
if stopInd(end) > size(data, 2)
    nEpochs = nEpochs-1;

end

% Initialize variables for outputs
nChans = size(data,1);  % number of EEG channels
SEF = nan(nChans,nEpochs);

% or define the variables one by one
[deltaDB,thetaDB,alphaDB,betaDB,broadDB] = deal(nan(nChans,nEpochs));

% Set 55 Hz cutoff for SEF calculation
cutoff = Find_Freq_Ind(55,epochLength);

% Set frequency ranges for EEG bands
deltaInd = Find_Freq_Ind([1 4],epochLength);
thetaInd = Find_Freq_Ind([4 8],epochLength);
alphaInd = Find_Freq_Ind([8 13],epochLength);
betaInd = Find_Freq_Ind([13 30],epochLength);




win = hann(windowLength);  % create hanning window based on number of time points in eegEpoch


% For each epoch, calculate the power spectrum
for chann = 1:nChans

    for epochInd = 1:nEpochs
        eegEpochtest(chann,epochInd,:) = data(chann,startInd(epochInd):stopInd(epochInd))'; 

        eegEpoch = data(chann,startInd(epochInd):stopInd(epochInd))'; % transpose data so each channel is in a column
        % Apply windowing function to epoch of EEG data
        eegEpochWin =win.*eegEpoch;  % apply window to EEG in each channel

        % Calculate FFT and EEG power
        fftVals = fft(eegEpochWin);

        % Calculate power; dimension is # frequencies x # channels
        powerVals = abs(fftVals).^2;
        powerVals = powerVals(1:length(powerVals)/2+1, :);  % remove negative frequencies (two-sided power to one-sided power)
        powerVals(2:end-1) = 2*powerVals(2:end-1); % doubling the power except for DC and Nyquist

        % Calculate SEF
        SEF(chann,epochInd) = Calc_SEF(powerVals,fs,cutoff);  % dimension is 1 x channels

        % Calculate power in each frequency band
        deltaDB(chann,epochInd) = Calc_Band_Specific_Total_Power(powerVals, deltaInd, fs, win);
        thetaDB(chann,epochInd) = Calc_Band_Specific_Total_Power(powerVals, thetaInd, fs, win);
        alphaDB(chann,epochInd) = Calc_Band_Specific_Total_Power(powerVals, alphaInd, fs, win);
        betaDB(chann,epochInd)  = Calc_Band_Specific_Total_Power(powerVals, betaInd, fs, win);
        broadDB(chann,epochInd) = Calc_Band_Specific_Total_Power(powerVals, 'all', fs, win);
    end

end

return
%% Functions
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
frequencies = linspace(0, fs/2, length(powerVals));

% Calculate the cumulative sum of the power spectrum
cumulativePower = cumsum(powerSpec);

% Normalize the cumulative power to create a cumulative distribution function (CDF)
cumulativePowerCDF = cumulativePower / cumulativePower(end);



% Find the index where the CDF reaches or exceeds the 95th percentile
[~,percentileIndex] = min(abs(cumulativePowerCDF-0.95));

% Find the frequency corresponding to this index
SEF = frequencies(percentileIndex);


return



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

% Alternative1 calculation
fInd = round(freq * epochLength)+1;
% fInd = round(freq * N / fs) + 1;
% Alternative2 calculation for readability:
% nyqF = fs/2;
% fInd = floor((freq * nyqF * epochLength) / nyqF)
% Alternative 3
% Find_Freq_Ind = @(freq, epochLength) floor(freq*epochLength)
return


function powerDecibel = Calc_Band_Specific_Total_Power(powerVals, freqInd, fs, win)
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


% Extract power values for the desired frequency band
if strcmp(freqInd, 'all')
    powerMat = powerVals(:, :);
else
    powerMat = powerVals(freqInd(1):freqInd(2), :);  % dimension is frequencies x channels
end

% Convert power to power spectral density
S = sum(win.^2);  % scaling factor based on hanning window function
psdEEG = powerMat /( fs * S);  % Divide by fs and scaling factor to obtain power spectral density (power/Hz)

% Calculate total power in the desired frequency band and convert to
% decibels
powerLinear = trapz(psdEEG,1);  % calculate the area;
powerDecibel = 10*log10(powerLinear);  % convert total power to decibels
return
