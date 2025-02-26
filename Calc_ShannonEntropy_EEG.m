function entropy = calc_Entropy_EEG(EEG_filt,aveOptimal)
% Calculates the Shannon entropy value for each channel 
% Inputs:   EEG_filt: Filtered clean EEG data
%           aveOptimal: Optimal bin count to use (default 390), test on own
%                       data!. Use the Freedman-Diaconis on signal x to calculate: 
%                       nBins (max(x)-min(signal))/(2*IQR(x)*length(x)^(-1/3))
% Outputs:  entropy: Shannon entropy value for each channel
%
% Previously called: "calculateEntropy.m"
% Code used in Smith et al. 2021
% Rachel J. Smith (Lopouratory 2019)
%
%
% V.1.0: Comment code, generalized code to fit generic EEG datsets 
%        Remove 19 channel label setup from function (DH)
% V.1.1: Uploaded onto github (DH)

    %data = EEG_filt(index,:);
    nChann = size(EEG_filt,1);
    entropy = nan(1,nChann);
    
    for i = 1:nChann
        
        [N,~] = hist(EEG_filt(i,:),aveOptimal);
        probs = N./sum(N);
        entVal = -1.*sum(probs.*log2(probs+eps));

        entropy(i) = entVal; % replaces switch case
        
        % Removed Switch case for labels below
%         chanLabel = labels{i};
%             switch chanLabel
%                 case 'Fp1'
%                     entropy(1) = entVal;
%                 case 'FP1'
%                     entropy(1) = entVal;
%                 case 'Fp2'
%                     entropy(2) = entVal;
%                 case 'FP2'
%                     entropy(2) = entVal;
%                 case 'F3'
%                     entropy(3) = entVal;
%                 case 'F4'
%                     entropy(4) = entVal;
%                 case 'C3'
%                     entropy(5) = entVal;
%                 case 'C4'
%                     entropy(6) = entVal;
%                 case 'P3'
%                     entropy(7) = entVal;
%                 case 'P4'
%                     entropy(8) = entVal;
%                 case 'O1'
%                     entropy(9) = entVal;
%                 case 'O2'
%                     entropy(10) = entVal;
%                 case 'F7'
%                     entropy(11) = entVal;
%                 case 'F8'
%                     entropy(12) = entVal;
%                 case 'T3'
%                     entropy(13) = entVal;
%                 case 'T4'
%                     entropy(14) = entVal;
%                 case 'T5'
%                     entropy(15) = entVal;
%                 case 'T6'
%                     entropy(16) = entVal;
%                 case 'Fz'
%                     entropy(17) = entVal;
%                 case 'FZ'
%                     entropy(17) = entVal;
%                 case 'Cz'
%                     entropy(18) = entVal;
%                 case 'CZ'
%                     entropy(18) = entVal;
%                 case 'Pz'
%                     entropy(19) = entVal;
%                 case 'PZ'
%                     entropy(19) = entVal;
%             end   
    end

end