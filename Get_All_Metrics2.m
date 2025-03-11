% By Venus 2.25.2025

clear all;close all;clc;

%select pre or post or both
dataDir = 'D:\Venus\Lab projects\PERF\UCB short clips';
phase = {'Pre','Post'};
status = {'Sleep1','Sleep2','Wake1','Wake2'};
headerName = 'reordered_hdr';
eegRecordName = 'reordered_record';
frequency = 'frequency';

desired_channel = 'Cz';

boxSize = 1;
order = 4;
delay = 1;



for p = 1: length(phase)

    for s = 1:length(status)

        currentDir = [dataDir '\' phase{p} '\' status{s}];
        count = 0;
        EEG_metric2 = {};
        fileNames = {};
        fileList = dir(fullfile(currentDir, '*.mat')); %reading the names of all .mat files

        for f = 1:length(fileList)
            patientMetrics = struct();
            
            fileNames{f,1} = fileList(f).name; % getting the name of the file we are analyzing
            
            %loading EEG signal and EEG header
            loadedData = load(fileList(f).name); % doesn't need to specify the path any more Load(currentDir ‘\’ fileNames{f}) 
            recordEEG = loadedData.(eegRecordName);
            hdrEEG = loadedData.(headerName);
            
            % finding the location of the desired channel
            [~,channel] = max(strcmp(hdrEEG.label,desired_channel));

            fs = unique(hdrEEG.(frequency));
            %Amplitude and Power Specral Density and SEF
            filterType = 'broadband';
            rereferenceMethod = 'EAR';
            artifactDetection = 'True';
            epochLength = 30; %% do this once and then use this for the rest of the code
            rerefEEG = Rereference_EEG(recordEEG, hdrEEG, rereferenceMethod);
            filteredEEG = Filter_EEG(rerefEEG, fs, filterType);
            epochStart = Find_Clean_Indices(filteredEEG,fs, epochLength)';
            epochStop = epochStart + epochLength * fs - 1;
            numLargeEpochs = length(epochStart);


            %% Amplitude
            
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
            [shanEntBeta,permEntBeta] = deal(nan(numLargeEpochs,1));
            for epochId = 1:numLargeEpochs
                shanEntBeta(epochId,1) = Calc_ShannonEntropy_EEG(filteredEEG(Cz,epochStart(epochId):epochStop(epochId)),fs,aveOptimal,epochLength);
                permEntBeta(epochId,1) = Calc_PermutationEntropy_EEG(filteredEEG(Cz,epochStart(epochId):epochStop(epochId)),fs,order,delay,epochLength);
            end

             patientMetrics.shanEntBeta = [shanEntBeta,epochStart,epochStop];
             patientMetrics.permEntBeta = [permEntBeta,epochStart,epochStop];
             

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

