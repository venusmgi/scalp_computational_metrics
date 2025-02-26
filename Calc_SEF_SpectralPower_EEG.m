function [freq95, powerMatFull] = calc_SEF_EEG(EEG_filt,fs,epochLength)
% Calculates the spectral edge frequency value for each channel    
% Inputs:   EEG_filt: Filtered clean EEG data
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

    nSecs = floor(size(EEG_filt,2)./fs);
    nEpochs = floor(nSecs./epochLength);
    nChan = size(EEG_filt,1);

    % Downsample to 200 hz
    if fs~=200
        EEG_new = nan(nChan,nSecs*200);
        for p = 1:nChan
            EEG_new(p,:) = interp1(1/fs:1/fs:nSecs,EEG_filt(p,1:nSecs*fs),1/200:1/200:nSecs);
        end
        EEG_filt = EEG_new;
        fs = 200;
    end
    
    % Initialize variables
    powerMat = nan(nChan,(epochLength*200)/2+1,nEpochs);   
    EEG = EEG_filt;
    
%     EEG =nan(size(EEG_filt));   
%     REMOVED: Add in NaNs for all the removed channels
%     for c = 1:size(labels,2)
%         chanLabel = labels{c};                
%         switch chanLabel
%             case 'Fp1'
%                 EEG(1,:) = EEG_filt(c,:);
%             case 'FP1'
%                 EEG(1,:) = EEG_filt(c,:);
%             case 'Fp2'
%                 EEG(2,:) = EEG_filt(c,:);
%             case 'FP2'
%                 EEG(2,:) = EEG_filt(c,:);
%             case 'F3'
%                 EEG(3,:) = EEG_filt(c,:);
%             case 'F4'
%                 EEG(4,:) = EEG_filt(c,:);
%             case 'C3'
%                 EEG(5,:) = EEG_filt(c,:);
%             case 'C4'
%                 EEG(6,:) = EEG_filt(c,:);
%             case 'P3'
%                 EEG(7,:) = EEG_filt(c,:);
%             case 'P4'
%                 EEG(8,:) = EEG_filt(c,:);
%             case 'O1'
%                 EEG(9,:) = EEG_filt(c,:);
%             case 'O2'
%                 EEG(10,:) = EEG_filt(c,:);
%             case 'F7'
%                 EEG(11,:) = EEG_filt(c,:);
%             case 'F8'
%                 EEG(12,:) = EEG_filt(c,:);
%             case 'T3'
%                 EEG(13,:) = EEG_filt(c,:);
%             case 'T4'
%                 EEG(14,:) = EEG_filt(c,:);
%             case 'T5'
%                 EEG(15,:) = EEG_filt(c,:);
%             case 'T6'
%                 EEG(16,:) = EEG_filt(c,:);
%             case 'Fz'
%                 EEG(17,:) = EEG_filt(c,:);
%             case 'FZ'
%                 EEG(17,:) = EEG_filt(c,:);
%             case 'Cz'
%                 EEG(18,:) = EEG_filt(c,:);
%             case 'CZ'
%                 EEG(18,:) = EEG_filt(c,:);
%             case 'Pz'
%                 EEG(19,:) = EEG_filt(c,:);
%             case 'PZ'
%                 EEG(19,:) = EEG_filt(c,:);
%             end   
%         end

    % For each epoch, calculate the power in each frequency point
    for n = 1:nEpochs
        patData = EEG(:,(n-1)*(epochLength*fs)+1:(n*epochLength*fs));
        for j = 1:size(EEG,1)
            channelData = patData(j,:);
            fftVals = abs(fft(channelData));
            powerMat(j,:,n) = fftVals(1:size(fftVals,2)/2+1);
        end
    end
    
    % Calculate sum, norm, and cumulative sum of powers
    meanPowerMatrix = (nanmean(powerMat,3)).^2;
    powerMatFull = meanPowerMatrix(:,1:275); % powerMatFull for SEF Full    
    meanPowerMatrix2 = meanPowerMatrix(:,1:275); % this is 55 Hz
    sumPower = sum(meanPowerMatrix2,2);
    normPower = meanPowerMatrix2./repmat(sumPower,1,275); %55 Hz
    cumSumNorm = cumsum(normPower,2);
    freq95 = nan(1,nChan);
    
    for k= 1:nChan
        [~,index] = min(abs(cumSumNorm(k,:)-0.95));
        fVec =  linspace(0,200/2,501);
        freq95(k) = fVec(index);
    end
    
    freq95(freq95==0) = nan;
        
end
