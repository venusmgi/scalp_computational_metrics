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
            num30SecEpochs = length(startingIndCleanEpoch);


            %% Amplitude
            amp = nan(num30SecEpochs,1);
            for epochId = 1:num30SecEpochs
                newEpochLength = 1;
                amp(epochId,1) = median(Calc_Amplitude_Range_EEG (filteredEEG(Cz,startingIndCleanEpoch(epochId):stopingIndCleanEpoch(epochId)), fs,newEpochLength));
            end
            
            amplitude_matrix = [amp,startingIndCleanEpoch,stopingIndCleanEpoch];
            patientMetrics.amplitude = amplitude_matrix;

            %% Power Specral Density and SEF
            [SEF,deltaDB,thetaDB,alphaDB,betaDB,broadDB] = deal(nan(num30SecEpochs,1));

            for epochId = 1:num30SecEpochs
                [tempSEF,tempDeltaDB,tempThetaDB,tempAlphaDB,tempBetaDB,tempBroadDB] = Calc_SEF_SpectralPower_EEG(filteredEEG(Cz,startingIndCleanEpoch(epochId):stopingIndCleanEpoch(epochId)),fs,newEpochLength);  

                SEF(epochId,:) = median(tempSEF);
                deltaDB(epochId,:) = median(tempDeltaDB);
                thetaDB(epochId,:) = median(tempThetaDB);
                alphaDB(epochId,:) = median(tempAlphaDB);
                betaDB(epochId,:) = median(tempBetaDB);
                broadDB(epochId,:) = median(tempBroadDB);

            end

            
            patientMetrics.psdSEF = [SEF,startingIndCleanEpoch,stopingIndCleanEpoch];
            patientMetrics.deltaPSD = [deltaDB,startingIndCleanEpoch,stopingIndCleanEpoch];
            patientMetrics.thetaPSD = [thetaDB,startingIndCleanEpoch,stopingIndCleanEpoch];
            patientMetrics.alphaPSD = [alphaDB,startingIndCleanEpoch,stopingIndCleanEpoch];
            patientMetrics.betaPSD = [betaDB,startingIndCleanEpoch,stopingIndCleanEpoch];
            patientMetrics.broadPSD = [broadDB,startingIndCleanEpoch,stopingIndCleanEpoch];

           %% Entropies 
            newEpochLength = 30;
            boxSize = 1;
            order = 4;
            delay = 1;

            % delta 
            filterType = 'delta';
            filteredEEG = Rereference_Filter_DetectArts(reordered_record,reordered_hdr,rereferenceMethod,filterType,epochLength,'False');
            aveOptimal =631;
            [shanEntDelta,permEntDelta] = deal(nan(num30SecEpochs,1));
            for epochId = 1:num30SecEpochs
                shanEntDelta(epochId,1) = Calc_ShannonEntropy_EEG(filteredEEG(Cz,startingIndCleanEpoch(epochId):stopingIndCleanEpoch(epochId)),fs,aveOptimal,epochLength);
                permEntDelta(epochId,1) = Calc_PermutationEntropy_EEG(filteredEEG(Cz,startingIndCleanEpoch(epochId):stopingIndCleanEpoch(epochId)),fs,order,delay,epochLength);
            end

            patientMetrics.shanEntDelta = [shanEntDelta,startingIndCleanEpoch,stopingIndCleanEpoch];
            patientMetrics.permEntDelta = [permEntDelta,startingIndCleanEpoch,stopingIndCleanEpoch];



            % theta 
            filterType = 'theta';
            filteredEEG = Rereference_Filter_DetectArts(reordered_record,reordered_hdr,rereferenceMethod,filterType,newEpochLength,'False');
            aveOptimal =735;
            [shanEntTheta,permEntTheta] = deal(nan(num30SecEpochs,1));
            for epochId = 1:num30SecEpochs
                shanEntTheta(epochId,1) = Calc_ShannonEntropy_EEG(filteredEEG(Cz,startingIndCleanEpoch(epochId):stopingIndCleanEpoch(epochId)),fs,aveOptimal,epochLength);
                permEntTheta(epochId,1) = Calc_PermutationEntropy_EEG(filteredEEG(Cz,startingIndCleanEpoch(epochId):stopingIndCleanEpoch(epochId)),fs,order,delay,epochLength);
            end

             patientMetrics.shanEntTheta = [shanEntTheta,startingIndCleanEpoch,stopingIndCleanEpoch];
             patientMetrics.permEntTheta = [permEntTheta,startingIndCleanEpoch,stopingIndCleanEpoch];

            % alpha 
            filterType = 'alpha';
            filteredEEG = Rereference_Filter_DetectArts(reordered_record,reordered_hdr,rereferenceMethod,filterType,newEpochLength,'False');
            aveOptimal =1008;
            [shanEntAlpha,permEntAlpha] = deal(nan(num30SecEpochs,1));
            for epochId = 1:num30SecEpochs
                shanEntAlpha(epochId,1) = Calc_ShannonEntropy_EEG(filteredEEG(Cz,startingIndCleanEpoch(epochId):stopingIndCleanEpoch(epochId)),fs,aveOptimal,epochLength);
                permEntAlpha(epochId,1) = Calc_PermutationEntropy_EEG(filteredEEG(Cz,startingIndCleanEpoch(epochId):stopingIndCleanEpoch(epochId)),fs,order,delay,epochLength);
            end

             patientMetrics.shanEntAlpha = [shanEntAlpha,startingIndCleanEpoch,stopingIndCleanEpoch];
             patientMetrics.permEntAlpha = [permEntAlpha,startingIndCleanEpoch,stopingIndCleanEpoch];
             

            % beta 
            filterType = 'beta';
            filteredEEG = Rereference_Filter_DetectArts(reordered_record,reordered_hdr,rereferenceMethod,filterType,newEpochLength,'False');
            aveOptimal =1517;
            [shanEntBeta,permEntBeta] = deal(nan(num30SecEpochs,1));
            for epochId = 1:num30SecEpochs
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

