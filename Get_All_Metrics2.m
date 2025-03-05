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
        EEG_metric2 = {};
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

            epochRatio = floor(epochLength/newEpochLength);
            amp = nan(30,length(startingIndCleanEpoch));
            for epochId = 1:length(startingIndCleanEpoch)
                newEpochLength = 1;
                
                amp(:,epochId) = Calc_Amplitude_Range_EEG (filteredEEG(Cz,startingIndCleanEpoch(epochId):stopingIndCleanEpoch(epochId)), fs,newEpochLength);
            end
            amp = reshape(amp,[],1);
            
            amplitude_matrix = [amp,newStartIdx,newStopIdx];
            patientMetrics.amplitude = amplitude_matrix;

            %% Power Specral Density and SEF

            newEpochLength = 5;
            newStartIdx = Get_Sub_Clean_Epochs(startingIndCleanEpoch, fs, epochLength, newEpochLength);
            newStopIdx = newStartIdx+fs*newEpochLength-1;

            epochRatio = floor(epochLength/newEpochLength);
            [SEF,deltaDB,thetaDB,alphaDB,betaDB,broadDB] = deal(nan(epochRatio,length(startingIndCleanEpoch)));
            
            for epochId = 1:length(startingIndCleanEpoch)
                [SEF(:,epochId),deltaDB(:,epochId),thetaDB(:,epochId),alphaDB(:,epochId),betaDB(:,epochId),broadDB(:,epochId)] = Calc_SEF_SpectralPower_EEG(filteredEEG(Cz,startingIndCleanEpoch(epochId):stopingIndCleanEpoch(epochId)),fs,newEpochLength);     
            end
            SEF = reshape(SEF,[],1);
            deltaDB = reshape(deltaDB,[],1);
            thetaDB = reshape(thetaDB,[],1);
            alphaDB = reshape(alphaDB,[],1);
            betaDB = reshape(betaDB,[],1);
            broadDB = reshape(broadDB,[],1);
            
            
            
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

            [shanEntDelta,permEntDelta] = deal(nan(length(startingIndCleanEpoch),1));
            for epochId = 1:length(startingIndCleanEpoch)
                shanEntDelta(epochId,1) = Calc_ShannonEntropy_EEG(filteredEEG(Cz,startingIndCleanEpoch(epochId):stopingIndCleanEpoch(epochId)),fs,aveOptimal,epochLength);
                permEntDelta(epochId,1) = Calc_PermutationEntropy_EEG(filteredEEG(Cz,startingIndCleanEpoch(epochId):stopingIndCleanEpoch(epochId)),fs,order,delay,epochLength);
            end
            patientMetrics.shanEntDelta = [shanEntDelta,startingIndCleanEpoch,stopingIndCleanEpoch];
            patientMetrics.permEntDelta = [permEntDelta,startingIndCleanEpoch,stopingIndCleanEpoch];



            % theta 
            filterType = 'theta';
            filteredEEG = Rereference_Filter_DetectArts(reordered_record,reordered_hdr,rereferenceMethod,filterType,newEpochLength,'False');
            aveOptimal =735;

            [shanEntTheta,permEntTheta] = deal(nan(length(startingIndCleanEpoch),1));
            for epochId = 1:length(startingIndCleanEpoch)
                shanEntTheta(epochId,1) = Calc_ShannonEntropy_EEG(filteredEEG(Cz,startingIndCleanEpoch(epochId):stopingIndCleanEpoch(epochId)),fs,aveOptimal,epochLength);
                permEntTheta(epochId,1) = Calc_PermutationEntropy_EEG(filteredEEG(Cz,startingIndCleanEpoch(epochId):stopingIndCleanEpoch(epochId)),fs,order,delay,epochLength);
            end
             patientMetrics.shanEntTheta = [shanEntTheta,startingIndCleanEpoch,stopingIndCleanEpoch];
             patientMetrics.permEntTheta = [permEntTheta,startingIndCleanEpoch,stopingIndCleanEpoch];

            % alpha 
            filterType = 'alpha';
            filteredEEG = Rereference_Filter_DetectArts(reordered_record,reordered_hdr,rereferenceMethod,filterType,newEpochLength,'False');
            aveOptimal =1008;

            [shanEntAlpha,permEntAlpha] = deal(nan(length(startingIndCleanEpoch),1));
            for epochId = 1:length(startingIndCleanEpoch)
                shanEntAlpha(epochId,1) = Calc_ShannonEntropy_EEG(filteredEEG(Cz,startingIndCleanEpoch(epochId):stopingIndCleanEpoch(epochId)),fs,aveOptimal,epochLength);
                permEntAlpha(epochId,1) = Calc_PermutationEntropy_EEG(filteredEEG(Cz,startingIndCleanEpoch(epochId):stopingIndCleanEpoch(epochId)),fs,order,delay,epochLength);
            end
             patientMetrics.shanEntAlpha = [shanEntAlpha,startingIndCleanEpoch,stopingIndCleanEpoch];
             patientMetrics.permEntAlpha = [permEntAlpha,startingIndCleanEpoch,stopingIndCleanEpoch];
             

            % beta 
            filterType = 'beta';
            filteredEEG = Rereference_Filter_DetectArts(reordered_record,reordered_hdr,rereferenceMethod,filterType,newEpochLength,'False');
            aveOptimal =1517;

            [shanEntBeta,permEntBeta] = deal(nan(length(startingIndCleanEpoch),1));
            for epochId = 1:length(startingIndCleanEpoch)
                shanEntBeta(epochId,1) = Calc_ShannonEntropy_EEG(filteredEEG(Cz,startingIndCleanEpoch(epochId):stopingIndCleanEpoch(epochId)),fs,aveOptimal,epochLength);
                permEntBeta(epochId,1) = Calc_PermutationEntropy_EEG(filteredEEG(Cz,startingIndCleanEpoch(epochId):stopingIndCleanEpoch(epochId)),fs,order,delay,epochLength);
            end
             patientMetrics.shanEntBeta = [shanEntBeta,startingIndCleanEpoch,stopingIndCleanEpoch];
             patientMetrics.permEntBeta = [permEntBeta,startingIndCleanEpoch,stopingIndCleanEpoch];


            EEG_metric2{count,1} = patientMetrics;
        end
        % Construct the filename using sprintf
        EEG_metric_filename = sprintf('%s_%s_results.mat', status{s}, phase{p});
        fileNames_filename = sprintf('%s_%s_filenames.mat', status{s}, phase{p});

        % Save the variables to the respective files
        save(EEG_metric_filename, 'EEG_metric2', '-v7.3');
        save(fileNames_filename, 'fileNames', '-v7.3');

    end
end



% currentDir = 'D:\Venus\Lab projects\Pre\Sleep1\Data\Age-matched Cases'
% fileList1 = dir(fullfile(currentDir, '*.mat'));
%
% currentDir = 'D:\Venus\Lab projects\Post\Sleep 1\Data\Age-matched Cases'
% fileList2 = dir(fullfile(currentDir, '*.mat'));