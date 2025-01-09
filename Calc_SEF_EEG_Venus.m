function [SEF,deltaDB,thetaDB,alphaDB,betaDB,broadDB] = Calc_SEF_EEG_Venus(data,fs,epochLength)
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

% Initialize variables for outputs
nSecs = floor(size(data,2)./fs);  % duration of data matrix, in seconds
nEpochs = floor(nSecs./epochLength);  % number of epochs in data
nChan = size(data,1);  % number of EEG channels


SEF = nan(nEpochs,nChan);

%% I can also define this as a struct:
bandPowers = struct('delta', deltaDB, 'theta', thetaDB, 'alpha', alphaDB, 'beta', betaDB, 'broad', broadDB);

% or define the variables one by one
deltaDB = nan(nEpochs,nChan);
thetaDB = nan(nEpochs,nChan);
alphaDB = nan(nEpochs,nChan);
betaDB = nan(nEpochs,nChan);
broadDB = nan(nEpochs,nChan);

% Set 55 Hz cutoff for SEF calculation
cutoff = Find_Freq_Ind(55,fs,epochLength);

% Set frequency ranges for EEG bands
deltaInd = Find_Freq_Ind([1 4],fs,epochLength);
thetaInd = Find_Freq_Ind([4 8],fs,epochLength);
alphaInd = Find_Freq_Ind([8 13],fs,epochLength);
betaInd = Find_Freq_Ind([13 30],fs,epochLength);




% For each epoch, calculate the power spectrum
for epochInd = 1:nEpochs
    %% Select epoch
    startInd = (epochInd-1)*(epochLength*fs)+1;
    stopInd = (epochInd*epochLength*fs);
    eegEpoch = data(:,startInd:stopInd)'; % transpose data so each channel is in a column
    
    % Apply windowing function to epoch of EEG data
    %% VENUS CHECK THIS, I HAVEN'T DEBUGGED IT
    win = hann(stopInd-startInd+1);  % create hanning window based on number of time points in eegEpoch
    eegEpochWin = repmat(win,1,nChan).*eegEpoch;  % apply window to EEG in each channel

    % Calculate FFT and EEG power
    fftVals = abs(fft(eegEpochWin)).^2;  % Calculate power; dimension is # frequencies x # channels
    fftVals = fftVals(1:size(fftVals,2)/2+1, :);  % remove negative frequencies

    % Calculate SEF
    SEF(epochInd,:) = Calc_SEF(fftVals,cutoff);  % dimension is 1 x channels

    % Calculate power in each frequency band
    deltaDB(epochInd,:) = Calc_Total_Power(fftVals, deltaInd, epochLength, fs, win);
    thetaDB(epochInd,:) = Calc_Total_Power(fftVals, thetaInd, epochLength, fs, win);
    alphaDB(epochInd,:) = Calc_Total_Power(fftVals, alphaInd, epochLength, fs, win);
    betaDB(epochInd,:) = Calc_Total_Power(fftVals, betaInd, epochLength, fs, win);
    broadDB(epochInd,:) = Calc_Total_Power(fftVals, 'all', epochLength, fs, win);



end

return


% Function to calculate SEF for one epoch
function SEF = Calc_SEF(fftVals, cutoff)
% CALC_SEF Calculates the Spectral Edge Frequency (SEF) for one epoch
%
% Inputs:
%   fftVals - Matrix of FFT values, frequencies x channels
%   cutoff  - Index corresponding to the maximum frequency to consider (e.g., 55 Hz)
%
% Output:
%   SEF     - Spectral Edge Frequency for each channel (95th percentile)

    powerSpec = fftVals(1:cutoff,:);  % frequencies (up to 55 Hz) x channels
    SEF = prctile(powerSpec,95,1); % calculate 95th percentile for each column
return

%% I can also define it as below:
% Calc_SEF = @(fftVals, cutoff) prctile(fftvals(1:cutoff,:), 95, 1);


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
    % Alternative calculation for readability:
    % nyqF = fs/2;
    % fInd = floor((freq * nyqF * epochLength) / nyqF)
return

%% I can also define it as below:
% Find_Freq_Ind = @(freq, epochLength) floor(freq*epochLength)

function powerDecibel = Calc_Total_Power(fftVals, freqInd, epochLength, fs, win)
% CALC_TOTAL_POWER Calculates the total power in a specified frequency band
%
% Inputs:
%   fftVals     - Matrix of FFT values, frequencies x channels
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
% Note that freqInd can be either a scalar or a vector

    % Extract power values for the desired frequency band
    if strcmp(freqInd, 'all')
        powerMat = fftVals(:, :);
    else
        powerMat = fftVals(freqInd(1):freqInd(2), :);  % dimension is frequencies x channels
    end

    % Convert power to power spectral density
    S = sum(win.^2);  % scaling factor based on windowing function
    psdEEG = (2/(fs*S))*powerMat;  % Divide by fs and multiply by scaling factor to obtain power spectral density (power/Hz)
    
    % Calculate total power in the desired frequency band and convert to
    % decibels
    powerIntegral = sum(psdEEG,1)*(1/epochLength);  % calculate the area; sum the power and multiply by the frequency increment
    powerDecibel = 10*log10(powerIntegral);  % convert total power to decibels
return
