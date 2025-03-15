function [z,lag, partial]=Calc_Max_Cross_Correlation(fs, maxLag, data1, data2, partialOn, ref)
% CALC_MAX_CROSS_CORRELATION Calculates the maximum cross-correlation between EEG signals.
%
% This function computes the cross-correlation for each 1-second window of EEG data,
% based on the methodology described by Kramer et al.
% This function was previously was called crossCorrFn.
%
% INPUTS:
%   fs        - Sampling frequency.
%   maxLag    - Maximum number of lags for cross-correlation (e.g., 200ms).
%   data1     - EEG data matrix 1 (channels x time).Usually 1 second of data.
%   data2     - EEG data matrix 2 (channels x time).Usually 1 second of data.
%   partialOn - Binary flag; 1 to include partial correlation calculation.
%   ref       - Common reference signal for partial cross-correlation.
%
% OUTPUTS:
%   z         - Matrix of z-values for maximum cross-correlation (nElec x nElec).
%   lag       - Matrix of lag times at maximum cross-correlation.
%   partial   - Binary matrix; 1 if value removed due to high partial correlation.
%
% NOTE: For typical cross-correlation calculations, data1 and data2 will be identical.
% This function calculates cross-correlation between all channel pairs.
% Data1 and data2 may differ for permutation resampling calculations.


n=fs;  % sampling rate


% Number of electrodes and datapoints
nElec = size(data1,1);
nTime = size(data1,2);

% Initialize matrices for results
z=zeros(size(data1,1),size(data1,1)); % z-value for max cross-correlation
lag=zeros(size(data1,1),size(data1,1)); % lag time at max cross-correlation
partial = zeros(nElec, nElec); % Partial correlation flag


% p=zeros(size(data1,1),size(data1,1)); % p-value (unused)

% Normalize data to have zero mean and unit variance
data1 = data1 - repmat(mean(data1,2),1,nTime); % subtract mean
data1 = data1./(repmat(std(data1,0,2),1,nTime));  % divide by standard dev.
data2 = data2 - repmat(mean(data2,2),1,nTime);
data2 = data2./(repmat(std(data2,0,2),1,nTime));


% Calculate autocorrelation
autoCorr1 = zeros(nElec, 2*fs+1);
autoCorr2 = zeros(nElec, 2*fs+1);
for j=1:nElec
    autoCorr1(j,:)=xcorr(data1(j,:),data1(j,:),fs,'biased');
    autoCorr2(j,:)=xcorr(data2(j,:),data2(j,:),fs,'biased');
end

% Binary matrix: 1 if value removed due to high partial correlation
partial = zeros(nElec);

% Calculate all cross-correlations
for j=1:nElec
    for k=(j+1):nElec
        % Calculate cross correlation
        x1=data1(j,:);
        x2=data2(k,:);
        [CC,lagVec]=xcorr(data1(j,:),data2(k,:),maxLag,'biased');
        
        % Find maximum cross correlation and associated lag time
        [maxCC,lagCC] = max(abs(CC));
        lag(j,k) = (lagVec(lagCC));
        lag(k,j) = lag(j,k);
        
        % Take Fisher transformation
        FCC = 0.5*log((1+maxCC)/(1-maxCC));
        
        % Calculate partial correlation to estimate influence of common reference
        if partialOn==1
            Rp2=partialcorr(data1(j,:)',data2(k,:)',-ref');
            R2=corrcoef(data1(j,:)',data2(k,:)');
            
            % Not significant if difference between partial and regular
            % is too large
            if (R2(1,2)-Rp2)>= 0.25
                lag(j,k)=0;
                lag(k,j)=lag(j,k);
                partial(j,k) = 1;
                partial(k,j) = 1;
            end
        end
        
        
        % Fisher transform for autocorrelations
        a1=0.5.*log((1+autoCorr1(j,:))./(1-autoCorr1(j,:)));
        a2=0.5.*log((1+autoCorr2(j,:))./(1-autoCorr2(j,:)));
        
        % Estimate the variance based on the autocorrelations
        auto_var=dot(a1,a2)/(n-abs(lag(j,k)));
        
        % Calculate z-value and p-value
        zValue = FCC/sqrt(auto_var);
        z(j,k)= zValue;
        z(k,j)=z(j,k);

        
    end
end
return