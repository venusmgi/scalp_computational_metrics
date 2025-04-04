% Script Overview
% ---------------
% This script processes EEG data files stored in a specified directory.
% It performs the following tasks:
%
% 1. Directory Setup:
%    - Defines the directory structure and parameters for different experimental phases (e.g., Pre/Post) and states (e.g., Sleep/Wake).
%
% 2. File Processing:
%    - Iterates over each .mat file in the directory, loading the EEG data and header information.
%
% 3. Channel Selection:
%    - Identifies the desired EEG channel (e.g., 'Cz') for analysis.
%
% 4. EEG Preprocessing:
%    - Re-references the EEG data using either the 'EAR' or 'CAR' method.
%    - Filters the EEG data according to specified frequency bands (broadband, delta, theta, alpha, beta).
%
% 5. Epoch Identification:
%    - Identifies clean epochs in the EEG data for further analysis.
%
% 6. Amplitude and Power Spectral Analysis:
%    - Calculates the amplitude and power spectral density.
%    - Computes the spectral edge frequency (SEF) for each epoch.
%
% 7. Entropy Calculations:
%    - Calculates Shannon and permutation entropy for different frequency bands.
%
% 8. Data Storage:
%    - Stores the calculated metrics in a structured format and saves them to .mat files.
%
% Functions Called by the Script:
% - Rereference_EEG: Re-references the EEG data using specified methods.
% - Filter_EEG: Applies filtering to the EEG data.
% - Find_Clean_Indices: Identifies clean epochs in the EEG data.
% - Calc_Amplitude_Range_EEG: Calculates the amplitude range of the EEG data.
% - Calc_SEF_SpectralPower_EEG: Calculates spectral power and spectral edge frequency.
% - Calc_ShannonEntropy_EEG: Calculates Shannon entropy of the EEG data.
% - Calc_PermutationEntropy_EEG: Calculates permutation entropy of the EEG data.
%
% Author: Venus
% Date: 2.25.2025

clear variables; close all; clc;  % Clear workspace, close all figures, and clear command window

%% Define directories and parameters for processing
dataDir = 'D:\Data\UCB Data';  % Directory containing EEG data
dataDir = 'D:\Venus\Lab projects\PERF\UCB short clips'; 
phase = {'Pre', 'Post'};  % Phases of the treatment
state = {'Sleep1', 'Sleep2', 'Wake1', 'Wake2'};  % Different states or conditions
headerName = 'reordered_hdr';  % Variable name for EEG header in .mat files
eegRecordName = 'reordered_record';  % Variable name for EEG data in .mat files
frequencyName = 'frequency';  % Field name for sampling frequency in header

% EEG channels to include in the analysis; write the list in a cell {}
channelsToAnalyze = {'Cz', 'Fz'};
% channelsToAnalyze = {'Cz'};

% Epoch length parameter: this script will calculate one value of each
% metric for each EEG epoch of this duration, in seconds
epochLength = 30;

% Sub-epoch length parameters: Overall, we will calculate one value of each
% metric for each epoch, but in some cases we want to subdivide that larger
% epoch into smaller chunks to do the analysis; for example, we may want to
% calculate the PSD in 5-second sub-epochs and take the mean value over all
% 6 sub-epochs in a larger 30-second epoch
subEpochLengthAmp = 1; % sub-epoch duration for amplitude, in seconds
subEpochLengthPSD = 5; % sub-epoch duration for PSD, in seconds

% Parameters for artifact detection
stdAbove = 7.5;        % Standard deviation threshold for artifact detection
buffer = 0.9;          % Buffer around detected artifacts (in seconds)
nArtChans = 1;         % Number of channels that must exceed threshold for artifact detection

% Parameters for amplitude, SEF and power spectral density calculation
rereferenceMethod = 'EAR';  % choose 'EAR' for linked ears or 'CAR' for common average reference

% Parameters for permutation entropy calculations
boxSize = 1;
order = 4;
delay = 1;

% Parameters for Shannon entropy, taken from Smith et al. 2021
% These values are calculated for a 30-second epoch length; if you would like
% to use a different epoch length to calculate Shannon entropy, you should
% calculate these values again using the Freedman-Diaconis rule:
%    nBins = (max(x) - min(x)) / (2 * IQR(x) * length(x)^(-1/3)).
nBinsDelta = 631;
nBinsTheta = 735;
nBinsAlpha = 1008;
nBinsBeta = 1517;

%% Loop through all EEG files and calculate all metrics for each epoch
% Loop over each phase, e.g., pre- or post-treatment
for p = 1:length(phase)
    % Loop over each state, e.g., wakefulness or sleep
    for s = 1:length(state)

        % Get file names of EEG datasets in current directory
        currentDir = [dataDir '\' phase{p} '\' state{s} '\Cases']; % Path to current data directory
       
        %% Venus changed:
        % we don't have pre treatment controls, or post treatment controls,
        % so I don't think adding Cases makes sense

        currentDir = [dataDir '\' phase{p} '\' state{s}];
        fileCell = struct2cell(dir(fullfile(currentDir, '*.mat')));  % List all .mat files in the directory
        fileNames = fileCell(1,:)'; % the first field in the "dir" structure is the file name
       
        % Initialize cell to store results; each element of the cell will
        % be results for one EEG file
        eegComputationalMetrics = cell(length(fileNames),1);

        % Loop over each file in the directory
        for f = 1:length(fileNames)
            % Initialize structure to store metrics for each patient
            patientMetrics = struct();

            % Load EEG signal and EEG header
            addpath(currentDir)
            loadedData = load(fileNames{f}); % doesn't need to specify the path any more Load(currentDir ‘\’ fileNames{f})
            recordEEG = loadedData.(eegRecordName);
            hdrEEG = loadedData.(headerName);

            % Find the index of each channel that will be analyzed
            nChan = length(channelsToAnalyze);
            chanVec = nan(1,nChan);
            for i = 1:nChan
                % This vector will contain a value of "1" at the index of the channel
                channelLoc = strcmp(hdrEEG.label,channelsToAnalyze{i});
                % If the channel is found, save the index; if not, throw error
                if sum(channelLoc)==0
                    error('Channel %s not found in EEG header.\n',channelsToAnalyze{i})
                else
                    chanVec(i) = find(channelLoc==1);
                end
            end
            
            % Extract the sampling rate and EEG duration
            fs = unique(hdrEEG.(frequencyName));
            N = size(recordEEG,2);

            % Re-reference EEG
            rerefEEG = Rereference_EEG(recordEEG, hdrEEG, rereferenceMethod); 


                 
            %% Why are we applying the filter after detecting the artifact? we decided to do it after filtering


            % Detect artifacts using the automated artifact detection function
            % NOTE: Input unfiltered EEG; function contains filtering
            % autoArts = get_automatedArtifacts_EEG(rerefEEG, fs, stdAbove, buffer, nArtChans);

             %%
            filteredEEG = Filter_EEG(rerefEEG, fs, 'broadband'); % Apply broadband filter to EEG
            autoArts = get_automatedArtifacts_EEG(filteredEEG, fs, stdAbove, buffer, nArtChans);

            epochStart = Find_Clean_Indices(N, fs, autoArts, epochLength); % Find start indices of the large clean epochs
            epochStop = epochStart + epochLength*fs - 1;  % Calculate stop indices for each of large clean epoch
            nEpoch = length(epochStart); % Number of large clean epochs
            
            % Save start and stop times of each epoch and list of channels
            patientMetrics.epochTimes = [epochStart,epochStop];
            patientMetrics.electrodes = channelsToAnalyze;
            
            %% Amplitude
            %% Why are we applying the filter after detecting the artifact? we decided to do it after filtering
            % % Apply broadband filter to EEG
            % filteredEEG = Filter_EEG(rerefEEG, fs, 'broadband'); % Apply broadband filter to EEG

            % Initialize vector for amplitude
            amp = nan(nEpoch,nChan);

            % Loop through each clean epoch
            for epochId = 1:nEpoch
                % Select the large epoch of EEG to analyze
                epochEEG = filteredEEG(chanVec,epochStart(epochId):epochStop(epochId));
                % Calculate amplitude; matrix will be nEpoch x nChan
                % Here we save the median across all sub-epochs
                amp(epochId,:) = median(Calc_Range_EEG(epochEEG, fs, subEpochLengthAmp),1);
            end

            patientMetrics.amplitude = amp;

            %% Power Spectral Density and SEF

            % Initialize matrices
            [SEF,deltaDB,thetaDB,alphaDB,betaDB,broadDB] = deal(nan(nEpoch,nChan));

            % Loop through all large epochs
            for epochId = 1:nEpoch
                % Calculate SEF and power for each sub-epoch
                epochEEG = filteredEEG(chanVec,epochStart(epochId):epochStop(epochId));
                [tempSEF,tempDeltaDB,tempThetaDB,tempAlphaDB,tempBetaDB,tempBroadDB]...
                    = Calc_SEF_SpectralPower_EEG(epochEEG,fs,subEpochLengthPSD);

                % Take the median value across sub-epochs
                SEF(epochId,:) = median(tempSEF,1);
                deltaDB(epochId,:) = median(tempDeltaDB,1);
                thetaDB(epochId,:) = median(tempThetaDB,1);
                alphaDB(epochId,:) = median(tempAlphaDB,1);
                betaDB(epochId,:) = median(tempBetaDB,1);
                broadDB(epochId,:) = median(tempBroadDB,1);
            end

            patientMetrics.SEF = SEF;
            patientMetrics.deltaPSD = deltaDB;
            patientMetrics.thetaPSD = thetaDB;
            patientMetrics.alphaPSD = alphaDB;
            patientMetrics.betaPSD = betaDB;
            patientMetrics.broadPSD = broadDB;

            %% Delta entropy
            % Delta frequency band
            filteredEEG = Filter_EEG(rerefEEG, fs, 'delta');

            % Calculate Shannon entropy and permutation entropy for each epoch
            entEEG = filteredEEG(chanVec,:);
            shanEntDelta = Calc_ShannonEntropy_EEG(entEEG,fs,nBinsDelta,epochLength,epochStart);
            permEntDelta = Calc_PermutationEntropy_EEG(entEEG,fs,order,delay,epochLength,epochStart);

            patientMetrics.shanEntDelta = shanEntDelta;
            patientMetrics.permEntDelta = permEntDelta;

            %% Theta entropy
            % Theta frequency band
            filteredEEG = Filter_EEG(rerefEEG, fs, 'theta');

            % Calculate Shannon entropy and permutation entropy for each epoch
            entEEG = filteredEEG(chanVec,:);
            shanEntTheta = Calc_ShannonEntropy_EEG(entEEG,fs,nBinsTheta,epochLength,epochStart);
            permEntTheta = Calc_PermutationEntropy_EEG(entEEG,fs,order,delay,epochLength,epochStart);

            patientMetrics.shanEntTheta = shanEntTheta;
            patientMetrics.permEntTheta = permEntTheta;

            %% Alpha entropy
            % Alpha frequency band
            filteredEEG = Filter_EEG(rerefEEG, fs, 'alpha');

            % Calculate Shannon entropy and permutation entropy for each epoch
            entEEG = filteredEEG(chanVec,:);
            shanEntAlpha = Calc_ShannonEntropy_EEG(entEEG,fs,nBinsAlpha,epochLength,epochStart);
            permEntAlpha = Calc_PermutationEntropy_EEG(entEEG,fs,order,delay,epochLength,epochStart);

            patientMetrics.shanEntAlpha = shanEntAlpha;
            patientMetrics.permEntAlpha = permEntAlpha;

            %% Beta entropy
            % Beta frequency band
            filteredEEG = Filter_EEG(rerefEEG, fs, 'beta');

            % Calculate Shannon entropy and permutation entropy for each epoch
            entEEG = filteredEEG(chanVec,:);
            shanEntBeta = Calc_ShannonEntropy_EEG(entEEG,fs,nBinsBeta,epochLength,epochStart);
            permEntBeta = Calc_PermutationEntropy_EEG(entEEG,fs,order,delay,epochLength,epochStart);
            
            patientMetrics.shanEntBeta = shanEntBeta;
            patientMetrics.permEntBeta = permEntBeta;

            % Store all metrics for each EEG file
            eegComputationalMetrics{f,1} = patientMetrics;

        end

        % Construct the filename using sprintf
        eegMetricFilename = sprintf('%s_%s_results.mat', state{s}, phase{p});

        % Save the variables to the respective files
        save(eegMetricFilename, "eegComputationalMetrics", "fileNames", '-v7.3');
    end
end
