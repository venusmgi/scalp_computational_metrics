% By Venus 2.25.2025

clear all;close all;clc;

%select pre or post or both
dataDir = 'D:\Venus\Lab projects\PERF\UCB short clips';
phase = {'Pre','Post'};
status = {'Sleep1','Sleep2','Wake1','Wake2'};


Cz = 18;

for p = 1: length(phase)
    phaseDir = [dataDir '\' phase{p}];
    for s = 1:length(status)
        currentDir = [phaseDir '\' status{s}];
        addpath(currentDir)
        count = 0;
        EEG_metric = {};
        fileNames = {};
        fileList = dir(fullfile(currentDir, '*.mat'));
        for f = 1:length(fileList)
            count = count+ 1;
            patientMetrics = struct();
            fileNames{count,1} = fileList(f).name; % getting the name of the file we are analyzing
            
            load(fileList(f).name)
            fs = reordered_hdr.frequency(1);
            %Amplitude and Power Specral Density and SEF
            filterType = 'broadband';
            rereferenceMethod = 'EAR';
            artifactDetection = 'True';
            epochLength = 30; %% do this once and then use this for the rest of the code


            [filteredEEG,startingIndCleanEpoch] = Rereference_Filter_DetectArts(reordered_record,reordered_hdr,rereferenceMethod,filterType,epochLength,artifactDetection);
            stopingIndCleanEpoch = startingIndCleanEpoch + epochLength * fs - 1;


            %% Amplitude
            newEpochLength = 1;
            newStartIdx = Get_Sub_Clean_Epochs(startingIndCleanEpoch, fs, epochLength, newEpochLength);
            newStopIdx = newStartIdx+fs*newEpochLength-1;

            % you can either use this method to calculate amplitude
            amp = Calc_Amplitude_Range_EEG (filteredEEG(Cz,:), fs,newEpochLength,newStartIdx)';
            
            amplitude_matrix = [amp,newStartIdx,newStopIdx];
            patientMetrics.amplitude = amplitude_matrix;

            %% Power Specral Density and SEF

            newEpochLength = 5;
            newStartIdx = Get_Sub_Clean_Epochs(startingIndCleanEpoch, fs, epochLength, newEpochLength);
            newStopIdx = newStartIdx+fs*newEpochLength-1;

            [SEF,deltaDB,thetaDB,alphaDB,betaDB,broadDB] = Calc_SEF_SpectralPower_EEG(filteredEEG(Cz,:),fs,newEpochLength,newStartIdx);
            
            patientMetrics.psdSEF = [SEF,newStartIdx,newStopIdx];
            patientMetrics.deltaPSD = [deltaDB,newStartIdx,newStopIdx];
            patientMetrics.thetaPSD = [thetaDB,newStartIdx,newStopIdx];
            patientMetrics.alphaPSD = [alphaDB,newStartIdx,newStopIdx];
            patientMetrics.betaPSD = [betaDB,newStartIdx,newStopIdx];
            patientMetrics.broadPSD = [broadDB,newStartIdx,newStopIdx];

           %% Entropies 
            newEpochLength = 30;
            boxSize = 1;
            order = 4;
            delay = 1;
            % delta 
            filterType = 'delta';
            filteredEEG = Rereference_Filter_DetectArts(reordered_record,reordered_hdr,rereferenceMethod,filterType,epochLength,'False');

            aveOptimal =631;
            shanEntDelta = Calc_ShannonEntropy_EEG(filteredEEG(Cz,:),fs,aveOptimal,epochLength,startingIndCleanEpoch);
            patientMetrics.shanEntDelta = [shanEntDelta,startingIndCleanEpoch,stopingIndCleanEpoch];

            permEntDelta = Calc_PermutationEntropy_EEG(filteredEEG(Cz,:),fs,order,delay,epochLength,startingIndCleanEpoch);
            patientMetrics.permEntDelta = [permEntDelta,startingIndCleanEpoch,stopingIndCleanEpoch];



            % theta 
            filterType = 'theta';
            filteredEEG = Rereference_Filter_DetectArts(reordered_record,reordered_hdr,rereferenceMethod,filterType,newEpochLength,'False');
            aveOptimal =735;
            shanEntTheta = Calc_ShannonEntropy_EEG(filteredEEG(Cz,:),fs,aveOptimal,epochLength,startingIndCleanEpoch);
            patientMetrics.shanEntTheta = [shanEntTheta,startingIndCleanEpoch,stopingIndCleanEpoch];

            permEntTheta = Calc_PermutationEntropy_EEG(filteredEEG(Cz,:),fs,order,delay,epochLength,startingIndCleanEpoch);
            patientMetrics.permEntTheta = [permEntTheta,startingIndCleanEpoch,stopingIndCleanEpoch];


            % alpha 
            filterType = 'alpha';
            filteredEEG = Rereference_Filter_DetectArts(reordered_record,reordered_hdr,rereferenceMethod,filterType,newEpochLength,'False');
            aveOptimal =1008;

            shanEntAlpha = Calc_ShannonEntropy_EEG(filteredEEG(Cz,:),fs,aveOptimal,newEpochLength,startingIndCleanEpoch);
            patientMetrics.shanEntAlpha = [shanEntAlpha,startingIndCleanEpoch,stopingIndCleanEpoch];

            permEntAlpha = Calc_PermutationEntropy_EEG(filteredEEG(Cz,:),fs,order,delay,newEpochLength,startingIndCleanEpoch);
            patientMetrics.permEntAlpha = [permEntAlpha,startingIndCleanEpoch,stopingIndCleanEpoch];


            % beta 
            filterType = 'beta';
            filteredEEG = Rereference_Filter_DetectArts(reordered_record,reordered_hdr,rereferenceMethod,filterType,newEpochLength,'False');
            aveOptimal =1517;

            shanEntBeta = Calc_ShannonEntropy_EEG(filteredEEG(Cz,:),fs,aveOptimal,newEpochLength,startingIndCleanEpoch);
            patientMetrics.shanEntBeta = [shanEntBeta,startingIndCleanEpoch,stopingIndCleanEpoch];

            permEntBeta = Calc_PermutationEntropy_EEG(filteredEEG(Cz,:),fs,order,delay,newEpochLength,startingIndCleanEpoch);
            patientMetrics.permEntBeta = [permEntBeta,startingIndCleanEpoch,stopingIndCleanEpoch];



            EEG_metric{count,1} = patientMetrics;
        end
        % Construct the filename using sprintf
        EEG_metric_filename = sprintf('%s_%s_results.mat', status{s}, phase{p});
        fileNames_filename = sprintf('%s_%s_filenames.mat', status{s}, phase{p});

        % Save the variables to the respective files
        save(EEG_metric_filename, 'EEG_metric', '-v7.3');
        save(fileNames_filename, 'fileNames', '-v7.3');

    end
end



% currentDir = 'D:\Venus\Lab projects\Pre\Sleep1\Data\Age-matched Cases'
% fileList1 = dir(fullfile(currentDir, '*.mat'));
%
% currentDir = 'D:\Venus\Lab projects\Post\Sleep 1\Data\Age-matched Cases'
% fileList2 = dir(fullfile(currentDir, '*.mat'));