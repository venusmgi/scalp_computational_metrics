function [MI, MVL, isCorrect] = Calc_PAC(amplitudeSignal, phaseSignal, fs, epochLength, startingIndices)
% CALC_PAC Calculate Phase-Amplitude Coupling (PAC) using MI and MVL methods
%
% This function computes phase-amplitude coupling between two filtered signals
% using two established methods: Tort's Modulation Index (MI) and Canolty's 
% Mean Vector Length (MVL). PAC quantifies the relationship between the phase 
% of low-frequency oscillations and the amplitude of high-frequency oscillations.
%
% INPUTS:
%   amplitudeSignal   - [nChannels × nSamples] matrix of high-frequency signal
%                       (e.g., gamma band 35-70 Hz) for amplitude extraction
%   phaseSignal       - [nChannels × nSamples] matrix of low-frequency signal
%                       (e.g., delta band 3-4 Hz) for phase extraction
%   fs                - Sampling frequency in Hz (e.g., 1000)
%   epochLength       - Duration of each epoch in seconds (e.g., 2)
%   startingIndices   - (Optional) Custom starting sample indices for epochs.
%                       If not provided, epochs are created from beginning of signal.
%
% OUTPUTS:
%   MI                - [nEpochs × nChannels] matrix of Modulation Index values
%                       (Tort et al., 2010). Higher values indicate stronger PAC.
%   MVL               - [nEpochs × nChannels] matrix of Mean Vector Length values
%                       (Canolty et al., 2006). Range: 0 (no coupling) to 1 (perfect coupling).
%   isCorrect         - Scalar flag indicating if MI calculations passed validation
%                       (1 = correct, 0 = potential calculation error)
%
% EXAMPLE USAGE:
%   % Calculate PAC between delta phase and gamma amplitude
%   [MI_vals, MVL_vals, valid] = Calc_PAC(gammaSignal, deltaSignal, 1000, 2);
%
% REFERENCES:
%   - Tort et al. (2010). Measuring phase-amplitude coupling between neuronal 
%     oscillations of different frequencies. J Neurophysiol, 104(2):1195-210.
%   - Canolty et al. (2006). High gamma power is phase-locked to theta 
%     oscillations in human neocortex. Science, 313(5793):1626-8.
%
% DEPENDENCIES:
%   - Calc_Epoch_Indices: Function to calculate epoch boundaries
%   - Tort_MI: Function to calculate Modulation Index
%   - Canolty_MVL: Function to calculate Mean Vector Length
%   - hilbert: MATLAB built-in function for Hilbert transform
%
% See also: Tort_MI, Canolty_MVL, Calc_Epoch_Indices

%% Input validation and dimension extraction
nChan = size(amplitudeSignal,1);  % number of EEG channels
nSamp = size(amplitudeSignal,2);  % number of samples (time)
numBins = 18;  % Number of phase bins for MI calculation

if nargin == 5
    [nEpochs, startInd, stopInd] = Calc_Epoch_Indices(nSamp, fs, epochLength, startingIndices);
else
    [nEpochs, startInd, stopInd] = Calc_Epoch_Indices(nSamp, fs, epochLength);
end

% % Apply filtering
% dataNotchFiltered = notchFilter(data, fs);
% phaseSignal = deltaFilter(dataNotchFiltered);      % Delta band (3-4 Hz) for phase
% amplitudeSignal = gammaFilter(dataNotchFiltered);  % Gamma band (35-70 Hz) for amplitude

%% Calculate PAC metrics for each channel

% Initialize output arrays
MI = nan(nEpochs, nChan);
MVL = nan(nEpochs, nChan);
correctnessCheck = nan(nEpochs, nChan);


% process each channel
for epochIdx = 1:nEpochs
    % Process each epoch
    for  channelIdx = 1:nChan
    
        % Define sample range for current epoch
        startSample = startInd(epochIdx);
        endSample = stopInd(epochIdx);
        
        % Extract epoch data
        epochPhaseSignal = phaseSignal(channelIdx, startSample:endSample);
        epochAmplitudeSignal = amplitudeSignal(channelIdx, startSample:endSample);
        
        % Apply Hilbert transform to extract phase and amplitude
        hilbertDelta = hilbert(epochPhaseSignal);
        phaseValues = angle(hilbertDelta) * 180 / pi;  % Convert to degrees
        
        hilbertGamma = hilbert(epochAmplitudeSignal);
        amplitudeValues = abs(hilbertGamma);
        
        % Calculate PAC metrics
        [MI(epochIdx, channelIdx), correctnessCheck(epochIdx, channelIdx)] = Tort_MI(amplitudeValues, phaseValues, numBins);
        
        [MVL(epochIdx, channelIdx)] = Canolty_MVL(amplitudeValues, phaseValues);
    end
   


end
 isCorrect = all(correctnessCheck(:) == 1);





end

