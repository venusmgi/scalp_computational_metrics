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

% Define directories and parameters for processing
dataDir = 'D:\Venus\Lab projects\PERF\UCB short clips';  % Directory containing EEG data
phase = {'Pre', 'Post'};  % Phases of the experiment
status = {'Sleep1', 'Sleep2', 'Wake1', 'Wake2'};  % Different states or conditions
headerName = 'reordered_hdr';  % Variable name for EEG header in .mat files
eegRecordName = 'reordered_record';  % Variable name for EEG data in .mat files
frequency = 'frequency';  % Field name for sampling frequency in header
desired_channel = 'Cz';  % Channel of interest for analysis

% Parameters for entropy calculations
boxSize = 1;
order = 4;
delay = 1;


% Loop over each phase
for p = 1: length(phase)
    % Loop over each status
    for s = 1:length(status)

        currentDir = [dataDir '\' phase{p} '\' status{s}]; % Path to current data directory
        EEG_metric2 = {}; % Initialize cell array to store metrics
        fileNames = {}; % Initialize cell array to store filenames
        fileList = dir(fullfile(currentDir, '*.mat'));  % List all .mat files in the directory

        % Loop over each file in the directory
        for f = 1:length(fileList)

            patientMetrics = struct(); % Initialize structure to store metrics for each patient
            fileNames{f,1} = fileList(f).name; % getting the name of the file we are analyzing

            %loading EEG signal and EEG header
            loadedData = load(fileList(f).name); % doesn't need to specify the path any more Load(currentDir ‘\’ fileNames{f})
            recordEEG = loadedData.(eegRecordName);
            hdrEEG = loadedData.(headerName);

            % finding the location of the desired channel
            [~,channel] = max(strcmp(hdrEEG.label,desired_channel));
            fs = unique(hdrEEG.(frequency));


            % Preprocess EEG data
            filterType = 'broadband';
            rereferenceMethod = 'EAR';
            artifactDetection = 'True';
            epochLength = 30; % Length of the largest clean epoch in seconds
            rerefEEG = Rereference_EEG(recordEEG, hdrEEG, rereferenceMethod); % Re-reference EEG
            filteredEEG = Filter_EEG(rerefEEG, fs, filterType); % Filter EEG
            epochStart = Find_Clean_Indices(filteredEEG,fs, epochLength)'; % Find start indices of the largest clean epochs
            epochStop = epochStart + epochLength * fs - 1;  % Calculate stop indices for each of largest clean epochs
            numLargeEpochs = length(epochStart); % Number of largest  clean epochs


            %% Amplitude

            % % Using the Alternative Method
            % % subEpochLength = 1;
            % % newStartIdx = Get_Sub_Clean_Epochs(epochStart, fs, epochLength, subEpochLength);
            % % newStopIdx = newStartIdx+fs*subEpochLength-1;
            % % amp = Calc_Amplitude_Range_EEG (filteredEEG(Cz,:), fs,subEpochLength,newStartIdx)';


            amp = nan(numLargeEpochs,1);
            for epochId = 1:numLargeEpochs
                amp(epochId,1) = median(Calc_Amplitude_Range_EEG (filteredEEG(Cz,epochStart(epochId):epochStop(epochId)), fs,subEpochLength));
            end


            patientMetrics.amplitude = [amp,epochStart,epochStop];

            %% Power Specral Density and SEF

            % % Using the Alternative Method
            % % subEpochLength = 5;
            % % newStartIdx = Get_Sub_Clean_Epochs(epochStart, fs, epochLength, subEpochLength);
            % % newStopIdx = newStartIdx+fs*subEpochLength-1;
            % % [SEF,DeltaDB,ThetaDB,AlphaDB,BetaDB,BroadDB] = Calc_SEF_SpectralPower_EEG(filteredEEG(Cz,:),fs,subEpochLength,newStartIdx);
            % % % SEF= median(reshape(SEF,6,[]),1);
            % % % deltaDB = median(reshape(DeltaDB,6,[]),1);
            % % % thetaDB = median(reshape(ThetaDB,6,[]),1);
            % % % alphaDB= median(reshape(AlphaDB,6,[]),1);
            % % % betaDB = median(reshape(BetaDB,6,[]),1);
            % % % broadDB = median(reshape(BroadDB,6,[]),1);


            % Defining a smaller epoch size, that can be summed up to be
            % the same size of the largest clean epoch size% Defining a
            % smaller epoch size, that can be summed up to be the same
            % size of the largest clean epoch size
            % For example, here we are calculating the power spectral
            % densities for epochs of 5 second, and each 6 second would
            % make a 30 secon clean epoch, which was defned abobe
            subEpochLength = 5;

            [SEF,deltaDB,thetaDB,alphaDB,betaDB,broadDB] = deal(nan(numLargeEpochs,1));
            for epochId = 1:numLargeEpochs
                [tempSEF,tempDeltaDB,tempThetaDB,tempAlphaDB,tempBetaDB,tempBroadDB] = Calc_SEF_SpectralPower_EEG(filteredEEG(Cz,epochStart(epochId):epochStop(epochId)),fs,subEpochLength);

                SEF(epochId,:) = median(tempSEF);
                deltaDB(epochId,:) = median(tempDeltaDB);
                thetaDB(epochId,:) = median(tempThetaDB);
                alphaDB(epochId,:) = median(tempAlphaDB);
                betaDB(epochId,:) = median(tempBetaDB);
                broadDB(epochId,:) = median(tempBroadDB);

            end


            patientMetrics.psdSEF = [SEF,epochStart,epochStop];
            patientMetrics.deltaPSD = [deltaDB,epochStart,epochStop];
            patientMetrics.thetaPSD = [thetaDB,epochStart,epochStop];
            patientMetrics.alphaPSD = [alphaDB,epochStart,epochStop];
            patientMetrics.betaPSD = [betaDB,epochStart,epochStop];
            patientMetrics.broadPSD = [broadDB,epochStart,epochStop];

            %% Entropies

            % delta
            filterType = 'delta';
            filteredEEG = Filter_EEG(rerefEEG, fs, filterType);
            aveOptimal =631;
            % % Using the Alternative Method
            % % shanEntDelta = Calc_ShannonEntropy_EEG(filteredEEG(Cz,:),fs,aveOptimal,epochLength,startingIndCleanEpoch);
            % % permEntDelta = Calc_PermutationEntropy_EEG(filteredEEG(Cz,:),fs,order,delay,epochLength,startingIndCleanEpoch);


            [shanEntDelta,permEntDelta] = deal(nan(numLargeEpochs,1));
            for epochId = 1:numLargeEpochs
                shanEntDelta(epochId,1) = Calc_ShannonEntropy_EEG(filteredEEG(Cz,epochStart(epochId):epochStop(epochId)),fs,aveOptimal,epochLength);
                permEntDelta(epochId,1) = Calc_PermutationEntropy_EEG(filteredEEG(Cz,epochStart(epochId):epochStop(epochId)),fs,order,delay,epochLength);
            end

            patientMetrics.shanEntDelta = [shanEntDelta,epochStart,epochStop];
            patientMetrics.permEntDelta = [permEntDelta,epochStart,epochStop];



            % theta
            filterType = 'theta';
            filteredEEG = Filter_EEG(rerefEEG, fs, filterType);
            aveOptimal =735;
            % % Using the Alternative Method
            % % shanEntTheta = Calc_ShannonEntropy_EEG(filteredEEG(Cz,:),fs,aveOptimal,epochLength,startingIndCleanEpoch);
            % % permEntTheta = Calc_PermutationEntropy_EEG(filteredEEG(Cz,:),fs,order,delay,epochLength,startingIndCleanEpoch);


            [shanEntTheta,permEntTheta] = deal(nan(numLargeEpochs,1));
            for epochId = 1:numLargeEpochs
                shanEntTheta(epochId,1) = Calc_ShannonEntropy_EEG(filteredEEG(Cz,epochStart(epochId):epochStop(epochId)),fs,aveOptimal,epochLength);
                permEntTheta(epochId,1) = Calc_PermutationEntropy_EEG(filteredEEG(Cz,epochStart(epochId):epochStop(epochId)),fs,order,delay,epochLength);
            end

            patientMetrics.shanEntTheta = [shanEntTheta,epochStart,epochStop];
            patientMetrics.permEntTheta = [permEntTheta,epochStart,epochStop];

            % alpha
            filterType = 'alpha';
            filteredEEG = Filter_EEG(rerefEEG, fs, filterType);
            aveOptimal =1008;

            % % Using the Alternative Method
            % % shanEntAlpha = Calc_ShannonEntropy_EEG(filteredEEG(Cz,:),fs,aveOptimal,newEpochLength,startingIndCleanEpoch);
            % % permEntAlpha = Calc_PermutationEntropy_EEG(filteredEEG(Cz,:),fs,order,delay,newEpochLength,startingIndCleanEpoch);


            [shanEntAlpha,permEntAlpha] = deal(nan(numLargeEpochs,1));
            for epochId = 1:numLargeEpochs
                shanEntAlpha(epochId,1) = Calc_ShannonEntropy_EEG(filteredEEG(Cz,epochStart(epochId):epochStop(epochId)),fs,aveOptimal,epochLength);
                permEntAlpha(epochId,1) = Calc_PermutationEntropy_EEG(filteredEEG(Cz,epochStart(epochId):epochStop(epochId)),fs,order,delay,epochLength);
            end


            patientMetrics.shanEntAlpha = [shanEntAlpha,epochStart,epochStop];
            patientMetrics.permEntAlpha = [permEntAlpha,epochStart,epochStop];


            % beta
            filterType = 'beta';
            filteredEEG = Filter_EEG(rerefEEG, fs, filterType);
            aveOptimal =1517;

            % % Using the Alternative Method
            % % shanEntBeta = Calc_ShannonEntropy_EEG(filteredEEG(Cz,:),fs,aveOptimal,newEpochLength,startingIndCleanEpoch);
            % % permEntBeta = Calc_PermutationEntropy_EEG(filteredEEG(Cz,:),fs,order,delay,newEpochLength,startingIndCleanEpoch);


            [shanEntBeta,permEntBeta] = deal(nan(numLargeEpochs,1));
            for epochId = 1:numLargeEpochs
                shanEntBeta(epochId,1) = Calc_ShannonEntropy_EEG(filteredEEG(Cz,epochStart(epochId):epochStop(epochId)),fs,aveOptimal,epochLength);
                permEntBeta(epochId,1) = Calc_PermutationEntropy_EEG(filteredEEG(Cz,epochStart(epochId):epochStop(epochId)),fs,order,delay,epochLength);
            end



            patientMetrics.shanEntBeta = [shanEntBeta,epochStart,epochStop];
            patientMetrics.permEntBeta = [permEntBeta,epochStart,epochStop];

            % Store metrics for each patient
            EEG_metric2{f,1} = patientMetrics;
        end
        % Construct the filename using sprintf
        EEG_metric_filename = sprintf('%s_%s_results.mat', status{s}, phase{p});
        fileNames_filename = sprintf('%s_%s_filenames.mat', status{s}, phase{p});

        % Save the variables to the respective files
        save(EEG_metric_filename, 'EEG_metric2', '-v7.3');
        save(fileNames_filename, 'fileNames', '-v7.3');

    end
end

