function [SEF,deltaDB,broadDB] = Calc_SEF_EEG_Venus_codeReview121824(data,fs,epochLength)
% Calculates the spectral edge frequency value for each channel
% Inputs:   EEG_filt: EEG data, channels x time
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
% V.1.0: Comment code, generalized code to fit generic EEG datsets
%        Remove 19 channel label setup from function
%        Combined calculateSEF_UCB and calculateSEF_UCB_full together
%        (both are identical but separately for some reason)
% V.1.1: Uploaded onto github (DH)
% V.1.2 : Venus added powerMat as the output, set fs to 200, removed
% resampling

nSecs = floor(size(data,2)./fs);  % duration of data matrix, in seconds
nEpochs = floor(nSecs./epochLength);  % number of epochs in data
nChan = size(data,1);  % number of EEG channels

% Initialize variables for outputs
SEF = nan(nEpochs,nChan);
deltaDB = nan(nEpochs,nChan);

% Set 55 Hz cutoff for SEF calculation
cutoff = Find_Freq_Ind(55,epochLength);

% Set frequency ranges for EEG bands
deltaInd = Find_Freq_Ind([1 4],epochLength);
% thetaInd = ...  % VENUS ADD THIS

% For each epoch, calculate the power spectrum
for epochInd = 1:nEpochs
    % Select epoch
    startInd = (epochInd-1)*(epochLength*fs)+1;
    stopInd = (epochInd*epochLength*fs);
    eegEpoch = data(:,startInd:stopInd)'; % transpose data so each channel is in a column
    
    % Apply windowing function to epoch of EEG data
    % VENUS CHECK THIS, I HAVEN'T DEBUGGED IT
    win = hann(stopInd-startInd+1);  % create hanning window based on number of time points in eegEpoch
    eegEpochWin = repmat(win,1,nChan).*eegEpoch;  % apply window to EEG in each channel

    % Calculate FFT and EEG power
    fftVals = abs(fft(eegEpochWin)).^2;  % Calculate power; dimension is # frequencies x # channels
    fftVals = fftVals(1:length(fftVals)/2+1, :);  % remove negative frequencies

    % Calculate SEF
    SEF(epochInd,:) = Calc_SEF(fftVals,cutoff);  % dimension is 1 x channels

    % Calculate power in each frequency band
    deltaDB(epochInd,:) = Calc_Total_Power(fftVals, deltaInd, epochLength, fs, win);
     broadDB(epochInd,:) = Calc_Total_Power(fftVals, 'all', epochLength, fs, win);

end

return


% Function to calculate SEF for one epoch
function SEF = Calc_SEF(fftVals, cutoff)
% ADD DESCRIPTIONS OF INPUTS AND OUTPUTS
    powerSpec = fftVals(1:cutoff,:);  % frequencies (up to 55 Hz) x channels
    SEF = prctile(powerSpec,95,1); % calculate 95th percentile for each column
return

% Function to find the index for a specific frequency value associated with the fft
% result
function fInd = Find_Freq_Ind(freq,epochLength)
% ADD DESCRIPTIONS OF INPUTS AND OUTPUTS
    fInd = floor(freq*epochLength);  %VENUS CHECK THIS!
return

function powerDecibel = Calc_Total_Power(fftVals, freqInd, epochLength, fs, win)
% ADD DESCRIPTIONS OF INPUTS and OUTPUTS
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


