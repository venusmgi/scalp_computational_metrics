% Script Overview
% ---------------
% This script processes EEG data files stored in a specified directory.
% It performs the following tasks:
%
% 1. Directory Setup:
%    - Defines the directory structure and parameters for different experimental phases and states (e.g., Pre, Post, Sleep, Wake).
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

clear all; close all; clc;  % Clear workspace, close all figures, and clear command window

%% Define directories and parameters for processing
dataDir = 'D:\Venus\Lab projects\PERF\UCB short clips';  % Directory containing EEG data
phase = {'Pre', 'Post'};  % Phases of the experiment
status = {'Sleep1', 'Sleep2', 'Wake1', 'Wake2'};  % Different states or conditions
headerName = 'reordered_hdr';  % Variable name for EEG header in .mat files
eegRecordName = 'reordered_record';  % Variable name for EEG data in .mat files
frequency = 'frequency';  % Field name for sampling frequency in header
desired_channel = {'Cz'};  % Channel of interest for analysis, make sure it is in a {}
 % number of channels we are analyzing


% patrameter for amplitude, SEF and power spectral denisty calculation
filterTypeBroad = 'broadband';
rereferenceMethod = 'EAR';
largestEpochLength = 30; % Length of the largest clean epoch in seconds

% Defining a smaller epoch size, that can be summed up to be
% the same size of the largest clean epoch size
% For example, here we are calculating the power spectral
% densities for epochs of 5 second, and each 6 second would
% make a 30 second clean epoch, which was defined above
subEpochLengthAmp = 1; % smaller epoch length that we want to calculate the amplitude for
subEpochLengthPSD = 5;

% Parameters for entropy calculations
boxSize = 1;
order = 4;
delay = 1;

filterTypeDelta = 'delta';
filterTypeTheta = 'theta';
filterTypeAlpha = 'alpha';
filterTypeBeta = 'beta';

% These valus are calculated for the 30 second epoch length, if the epoch 
% length for which you want to calculate the shannon enptropy changes from 
% 30 second, you should calculate this value again using Freedman-Diaconis 
% rule on the signal to calculate:
% nBins = (max(x) - min(x)) / (2 * IQR(x) * length(x)^(-1/3)).
numBinsDelta = 631;
numBinsTheta = 735;
numBinsAlpha = 1008;
numBinsBeta = 1517;


% Loop over each phase
for p = 1: length(phase)
    % Loop over each status
    for s = 1:length(status)

        currentDir = [dataDir '\' phase{p} '\' status{s}]; % Path to current data directory
        eegComputationalMetrics = {}; % Initialize cell array to store metrics
        fileNames = {}; % Initialize cell array to store filenames
        fileList = dir(fullfile(currentDir, '*.mat'));  % List all .mat files in the directory

        % Loop over each file in the directory
        for f = 1:length(fileList)

            patientMetrics = struct(); % Initialize structure to store metrics for each patient
            fileNames{f,1} = fileList(f).name; % getting the name of the file we are analyzing

            %loading EEG signal and EEG header
            addpath(currentDir)
            loadedData = load(fileList(f).name); % doesn't need to specify the path any more Load(currentDir ‘\’ fileNames{f})
            recordEEG = loadedData.(eegRecordName);
            hdrEEG = loadedData.(headerName);

            % finding the location of the desired channel
            numChanns = length(desired_channel);
            channel = nan(1,numChanns);
            for i = 1:numChanns
                [~,channel(i)] = max(strcmp(hdrEEG.label,desired_channel{i}));
            end
            fs = unique(hdrEEG.(frequency));


            % Preprocess EEG data
            rerefEEG = Rereference_EEG(recordEEG, hdrEEG, rereferenceMethod); % Re-reference EEG
            filteredEEG = Filter_EEG(rerefEEG, fs, filterTypeBroad); % Filter EEG
            epochStart = Find_Clean_Indices(filteredEEG,fs, largestEpochLength)'; % Find start indices of the largest clean epochs
            epochStop = epochStart + largestEpochLength * fs - 1;  % Calculate stop indices for each of largest clean epochs
            numLargeEpochs = length(epochStart); % Number of largest  clean epochs
            
            %% Amplitude

            amp = nan(numLargeEpochs,numChanns);
            for epochId = 1:numLargeEpochs
                %% Change the Calc_Amplitude_Range_EEG function name so user knows it is epoch based
                amp(epochId,:) = median(Calc_Amplitude_Range_EEG (filteredEEG(channel,epochStart(epochId):epochStop(epochId)), fs,subEpochLengthAmp),2);
            end

            patientMetrics.amplitude = [amp,epochStart,epochStop]; %what if we are analyzing for more than 1 channel?

            %% Power Specral Density and SEF

            [SEF,deltaDB,thetaDB,alphaDB,betaDB,broadDB] = deal(nan(numLargeEpochs,numChanns));
            for epochId = 1:numLargeEpochs
                [tempSEF,tempDeltaDB,tempThetaDB,tempAlphaDB,tempBetaDB,tempBroadDB] = Calc_SEF_SpectralPower_EEG(filteredEEG(channel,epochStart(epochId):epochStop(epochId)),fs,subEpochLengthPSD);

                SEF(epochId,:) = median(tempSEF,2); %median across epochs
                deltaDB(epochId,:) = median(tempDeltaDB,2);
                thetaDB(epochId,:) = median(tempThetaDB,2);
                alphaDB(epochId,:) = median(tempAlphaDB,2);
                betaDB(epochId,:) = median(tempBetaDB,2);
                broadDB(epochId,:) = median(tempBroadDB,2);

            end



            patientMetrics.SEF = [SEF,epochStart,epochStop];
            patientMetrics.deltaPSD = [deltaDB,epochStart,epochStop];
            patientMetrics.thetaPSD = [thetaDB,epochStart,epochStop];
            patientMetrics.alphaPSD = [alphaDB,epochStart,epochStop];
            patientMetrics.betaPSD = [betaDB,epochStart,epochStop];
            patientMetrics.broadPSD = [broadDB,epochStart,epochStop];

            %% Entropies

            % delta
            filteredEEG = Filter_EEG(rerefEEG, fs, filterTypeDelta);

            % % Using the Alternative Method
            % % shanEntDelta = Calc_ShannonEntropy_EEG(filteredEEG(channel,:),fs,numBinsDelta,largestEpochLength,epochStart);
            % % permEntDelta = Calc_PermutationEntropy_EEG(filteredEEG(channel,:),fs,order,delay,largestEpochLength,epochStart);
            [shanEntDelta,permEntDelta] = deal(nan(numLargeEpochs,numChanns));
            for epochId = 1:numLargeEpochs
                shanEntDelta(epochId,:) = Calc_ShannonEntropy_EEG(filteredEEG(channel,epochStart(epochId):epochStop(epochId)),fs,numBinsDelta,largestEpochLength);
                permEntDelta(epochId,:) = Calc_PermutationEntropy_EEG(filteredEEG(channel,epochStart(epochId):epochStop(epochId)),fs,order,delay,largestEpochLength);
            end

            patientMetrics.shanEntDelta = [shanEntDelta,epochStart,epochStop];
            patientMetrics.permEntDelta = [permEntDelta,epochStart,epochStop];

            % theta
            filteredEEG = Filter_EEG(rerefEEG, fs, filterTypeTheta);
            % % Using the Alternative Method
            % % shanEntTheta = Calc_ShannonEntropy_EEG(filteredEEG(channel,:),fs,numBinsTheta,largestEpochLength,epochStart);
            % % permEntTheta = Calc_PermutationEntropy_EEG(filteredEEG(channel,:),fs,order,delay,largestEpochLength,epochStart);


            [shanEntTheta,permEntTheta] = deal(nan(numLargeEpochs,numChanns));
            for epochId = 1:numLargeEpochs
                shanEntTheta(epochId,:) = Calc_ShannonEntropy_EEG(filteredEEG(channel,epochStart(epochId):epochStop(epochId)),fs,numBinsTheta,largestEpochLength);
                permEntTheta(epochId,:) = Calc_PermutationEntropy_EEG(filteredEEG(channel,epochStart(epochId):epochStop(epochId)),fs,order,delay,largestEpochLength);
            end

            patientMetrics.shanEntTheta = [shanEntTheta,epochStart,epochStop];
            patientMetrics.permEntTheta = [permEntTheta,epochStart,epochStop];

            % alpha
            filteredEEG = Filter_EEG(rerefEEG, fs, filterTypeAlpha);
            % % Using the Alternative Method
            % shanEntAlpha = Calc_ShannonEntropy_EEG(filteredEEG(channel,:),fs,numBinsAlpha,largestEpochLength,epochStart);
            % permEntAlpha = Calc_PermutationEntropy_EEG(filteredEEG(channel,:),fs,order,delay,largestEpochLength,epochStart);
            [shanEntAlpha,permEntAlpha] = deal(nan(numLargeEpochs,numChanns));
            for epochId = 1:numLargeEpochs
                shanEntAlpha(epochId,:) = Calc_ShannonEntropy_EEG(filteredEEG(channel,epochStart(epochId):epochStop(epochId)),fs,numBinsAlpha,largestEpochLength);
                permEntAlpha(epochId,:) = Calc_PermutationEntropy_EEG(filteredEEG(channel,epochStart(epochId):epochStop(epochId)),fs,order,delay,largestEpochLength);
            end


            patientMetrics.shanEntAlpha = [shanEntAlpha,epochStart,epochStop];
            patientMetrics.permEntAlpha = [permEntAlpha,epochStart,epochStop];


            % beta
            filteredEEG = Filter_EEG(rerefEEG, fs, filterTypeBeta);
            % % Using the Alternative Method
            % % shanEntBeta = Calc_ShannonEntropy_EEG(filteredEEG(channel,:),fs,numBinsBeta,largestEpochLength,epochStart);
            % % permEntBeta = Calc_PermutationEntropy_EEG(filteredEEG(channel,:),fs,order,delay,largestEpochLength,epochStart);
            [shanEntBeta,permEntBeta] = deal(nan(numLargeEpochs,numChanns));
            for epochId = 1:numLargeEpochs
                shanEntBeta(epochId,:) = Calc_ShannonEntropy_EEG(filteredEEG(channel,epochStart(epochId):epochStop(epochId)),fs,numBinsBeta,largestEpochLength);
                permEntBeta(epochId,:) = Calc_PermutationEntropy_EEG(filteredEEG(channel,epochStart(epochId):epochStop(epochId)),fs,order,delay,largestEpochLength);
            end

            patientMetrics.shanEntBeta = [shanEntBeta,epochStart,epochStop];
            patientMetrics.permEntBeta = [permEntBeta,epochStart,epochStop];

            % Store metrics for each patient
            eegComputationalMetrics{f,1} = patientMetrics;




        end
        % Construct the filename using sprintf
        eegMetricFilename = sprintf('%s_%s_results.mat', status{s}, phase{p});
        edfNames_filename = sprintf('%s_%s_filenames.mat', status{s}, phase{p});

        % Save the variables to the respective files
        save(eegMetricFilename, 'eegComputationalMetrics', '-v7.3');
        save(edfNames_filename, 'fileNames', '-v7.3');

    end
end

