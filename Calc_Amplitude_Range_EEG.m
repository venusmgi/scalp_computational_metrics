function [ampMatrix] = Calc_Amplitude_Range_EEG(data, fs)
% Outputs a matrix of all amplitude values for each one second window,
% calculated using the range (peak-to-peak).
% If artifactual epochs need to be removed, they can be identified using
% get_cleanData_EEG epochs and the artifactual amplitude values can be
% removed later.
%
% Inputs:
%    data: EEG data, channels x time points
%    fs: Sampling rate
% Outputs:
%    ampMatrix: matrix of all amplitude values for each one second window
%               (channels x nSeconds)
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
% - let user input length of window (do not restrict to 1 second); use
% variable "epochLength" to be consistent with SEF/Power spec code
% - transpose data matrix in main code for faster processing?

nChan = size(data,1);  % number of EEG channels
nSamps = size(data,2);  % number of samples (time)
nSecs = floor(nSamps./fs);  % number of 1-second windows

% Check that the data matrix is in the correct orientation
assert(nChan < nSamps, 'Warning: Number of channels is greater than number of time samples! Data matrix may need to be transposed.')

% Check that the data contains at least 1 second of data
assert(nSamps >= fs, 'Warning: Data matrix contains less than one second of data.')

% Initialize matrix for amplitude values, channels x 1-sec windows
ampMatrix = nan(nChan,nSecs);

% Loop through all 1-second windows
for win=1:nSecs
    startInd = (win-1)*fs+1;  % starting index for 1-second window
    stopInd = win*fs;  % ending index for 1-second window
    eegWin = data(:,startInd:stopInd);  % select 1-second window of EEG
    ampMatrix(:,win) = max(eegWin,[],2) - min(eegWin,[],2); % Calculate range for this window of data
end

