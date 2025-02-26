function [freq95, powerMatFull,powerMat] = calc_SEF_EEG_Venus(EEG_filt,fs,epochLength)
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
% V.1.2 : Venus added powerMat as the output, set fs to 200, removed
% resampling 


    nSecs = floor(size(EEG_filt,2)./fs);
    nEpochs = floor(nSecs./epochLength);
    nChan = size(EEG_filt,1);

    % % Downsample to 200 hz
    % if fs~=200
    %     EEG_new = nan(nChan,nSecs*20);
    %     for p = 1:nChan
    %         EEG_new(p,:) = interp1(1/fs:1/fs:nSecs,EEG_filt(p,1:nSecs*fs),1/20:1/20:nSecs);
    %     end
    %     EEG_filt = EEG_new;
    %     fs = 200;
    % end
    
    % Initialize variables
    powerMat = nan(nChan,(epochLength*200)/2+1,nEpochs);   
    EEG = EEG_filt;
    

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

    meanPowerMatrix = (mean(powerMat,3,'omitnan')).^2;
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
