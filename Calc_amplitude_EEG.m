function [ampMatrix] = Calc_amplitude_EEG(data,fs)
% Outputs a matrix of all amplitude values for each one second window
% The data should already be divided into one second windows from the
% get_cleanData_EEG epochs. 
% Inputs:   data: Filtered clean EEG data
%           fs: Sampling rate
% Outputs:  ampMatrix:  matrix of all amplitude values for each one second window
% (channels x nSeconds)
%
% Previously called: "calculate_Amplitude_full.m"
% Code used in Smith et al. 2021
% Rachel J. Smith (Lopouratory 2019)
%
% V.1.0: Comment code, generalized code to fit generic EEG datsets (DH)
%        Remove 19 channel label setup from function  
% V.1.1: Uploaded onto github (DH)


 nsamps = size(data,2);
 nSecs = floor(nsamps./fs);
 
   
   ampMatrix = nan(size(data,1),nSecs);
   
   for ch = 1:size(data,1)
       
        amps = zeros(1,nSecs);
        for j=1:nSecs
            EEG_vec = data(ch,(j-1)*fs+1:j*fs);
            amps(j) = max(EEG_vec)-min(EEG_vec);
        end
        
        ampMatrix(ch,:) = amps; % replaces switch case below
%        chanLabel = labels{ch};
%             switch chanLabel
%                 case 'Fp1'
%                     ampMatrix(1,:) = amps;
%                 case 'FP1'
%                     ampMatrix(1,:) = amps;
%                 case 'Fp2'
%                     ampMatrix(2,:) = amps;
%                 case 'FP2'
%                     ampMatrix(2,:) = amps;
%                 case 'F3'
%                     ampMatrix(3,:) = amps;
%                 case 'F4'
%                     ampMatrix(4,:) = amps;
%                 case 'C3'
%                     ampMatrix(5,:) = amps;
%                 case 'C4'
%                     ampMatrix(6,:) = amps;
%                 case 'P3'
%                     ampMatrix(7,:) = amps;
%                 case 'P4'
%                     ampMatrix(8,:) = amps;
%                 case 'O1'
%                     ampMatrix(9,:) = amps;
%                 case 'O2'
%                     ampMatrix(10,:) = amps;
%                 case 'F7'
%                     ampMatrix(11,:) = amps;
%                 case 'F8'
%                     ampMatrix(12,:) = amps;
%                 case 'T3'
%                     ampMatrix(13,:) = amps;
%                 case 'T4'
%                     ampMatrix(14,:) = amps;
%                 case 'T5'
%                     ampMatrix(15,:) = amps;
%                 case 'T6'
%                     ampMatrix(16,:) = amps;
%                 case 'Fz'
%                     ampMatrix(17,:) = amps;
%                 case 'FZ'
%                     ampMatrix(17,:) = amps;
%                 case 'Cz'
%                     ampMatrix(18,:) = amps;
%                 case 'CZ'
%                     ampMatrix(18,:) = amps;
%                 case 'Pz'
%                     ampMatrix(19,:) = amps;
%                 case 'PZ'
%                     ampMatrix(19,:) = amps;
%             end   
    end
   
end