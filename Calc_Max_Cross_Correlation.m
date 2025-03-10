% Calculate cross-correlation for each 1-second window of EEG; based on
% Kramer et al.
% Previously was called crossCorrFn
% INPUTS:
%     fs = sampling frequency
%     maxLag = maximum number of lags for cross-correlation (generally 200ms)
%     data1 = 1 second of data (for EEG, rows = channels, columns = time)
%     data2 = 1 second of data (for EEG, rows = channels, columns = time)
%     partialOn = binary; 1 will include partical correlation calculation
%     ref = common reference signal, for partial cross corelation
%
% OUTPUTS:
%     z = z-value for max cross-correlation; matrix (nElec x nElec)
%     lag = lag time at max cross-correlation
%     partial = binary matrix: 1 if value removed due to high partial
%               correlation (large influence of reference)
%
% NOTE: For a typical cross-correlation calculation, data1 and data2 will
% be exactly the same.  This function calculates cross correlation between
% all channel pairs. The matrices data1 and data2 will be different for the
% permutation resampling calculation, where one of the inputs represents
% time-shifted data.

function [z,lag, partial]=Calc_Max_Cross_Correlation(fs, maxLag, data1, data2, partialOn, ref)

n=fs;  % sampling rate
% n2=81; % this only shows up in calc for a and b (unused)

nElec = size(data1,1);
nTime = size(data1,2);

% Initialize matrices for results
z=zeros(size(data1,1),size(data1,1)); % z-value for max cross-correlation
lag=zeros(size(data1,1),size(data1,1)); % lag time at max cross-correlation
% p=zeros(size(data1,1),size(data1,1)); % p-value (unused)

% Ensure that the data has zero mean and unit variance
data1 = data1 - repmat(mean(data1,2),1,nTime); % subtract mean
data1 = data1./(repmat(std(data1,0,2),1,nTime));  % divide by standard dev.
data2 = data2 - repmat(mean(data2,2),1,nTime);
data2 = data2./(repmat(std(data2,0,2),1,nTime));

% Calculate autocorrelation
autoCorr1 = zeros(nElec, 2*fs+1);
autoCorr2 = zeros(nElec, 2*fs+1);
for j=1:nElec
    x1=data1(j,:);
    x2=data2(j,:);
    autoCorr1(j,:)=xcorr(x1,x1,fs,'biased');
    autoCorr2(j,:)=xcorr(x2,x2,fs,'biased');
end

% Binary matrix: 1 if value removed due to high partial correlation
partial = zeros(nElec);

% Calculate all cross-correlations
for j=1:nElec
    for k=(j+1):nElec
        % Calculate cross correlation
        x1=data1(j,:);
        x2=data2(k,:);
        [CC,lagVec]=xcorr(x1,x2,maxLag,'biased');
        
        % Find maximum cross correlation and associated lag time
        [maxCC,lagCC] = max(abs(CC));
        lag(j,k) = (lagVec(lagCC));
        lag(k,j) = lag(j,k);
        
        % Take Fisher transformation
        FCC = 0.5*log((1+maxCC)/(1-maxCC));
        
        % Calculate partial correlation to estimate influence of common reference
        if partialOn==1
            Rp2=partialcorr(x1',x2',-ref');
            R2=corrcoef(x1',x2');
            
            % Not significant if difference between partial and regular
            % is too large
            if (R2(1,2)-Rp2)>= 0.25
                lag(j,k)=0;
                lag(k,j)=lag(j,k);
                partial(j,k) = 1;
                partial(k,j) = 1;
            end
        end
        
        % STATISTICS -- from Kramer paper (unused)
        % Coefficients for the p-value calculation (we use permutation resampling instead)
        %             a=sqrt(2*log(n2));
        %             b=a-(log(log(n2))+log(4*pi))/(2*a);
        
        % Fisher transform for autocorrelations
        a1=0.5.*log((1+autoCorr1(j,:))./(1-autoCorr1(j,:)));
        a2=0.5.*log((1+autoCorr2(j,:))./(1-autoCorr2(j,:)));
        
        % Estimate the variance based on the autocorrelations
        auto_var=dot(a1,a2)/(n-abs(lag(j,k)));
        
        % Calculate z-value and p-value
        zValue = FCC/sqrt(auto_var);
        z(j,k)= zValue;
        z(k,j)=z(j,k);
        % From Kramer paper (unused; we did permutation resampling instead)
        % p(j,k)=(1-exp(-2*exp(-a*(zValue-b))));
        % p(k,j)=p(j,k);
        
    end
end
return