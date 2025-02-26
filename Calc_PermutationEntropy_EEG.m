function [entBroad,newSize] = calc_PermEntropy_EEG(EEG_filt,electrode,order,delay)
% Calculates the Permutation Entropy of the electrode chosen
% Inputs:   EEG_filt:   filtered clean time series;
%           electrode:  Electrode number of choice (default: CZ electrode)
%           order:      order of permuation entropy (default 5)
%           delay:      delay time of permuation entropy (default 1)
% Outputs:  entBroad:   Permutation entropy 
%           newSize:    Sample points of EEG_filt
%
% Previously called: "calculate_PermEntropy_Derek.m"
% Revised from code used in Smith et al. 2021
% Uses the pec function (Ouyang et al.)
% Derek Hu (Lopouratory, 2021)
%
% V.1.0: Comment code, generalized code to fit generic EEG datsets 
%        Removed the use of channel labels from function
%        Allow the choice of electrode to calculate pec of (CZ pref) (DH)
% V.1.1: Uploaded onto github (DH)

% findCz = strcmp(labels,'CZ');
% if sum(double(findCz))~=1
%     disp('THIS ONE DIDNT HAVE CZ');
% else
%     EEG_Cz = EEG_filt(findCz,:);

    EEG_elec = EEG_filt(electrode,:);
    [permEnt,~] = pec(EEG_elec,order,delay);
    entBroad = permEnt;
    
% end

    newSize = size(EEG_filt,2);

end
