function [EEG_reref] = Rereference_EEG(EEG_record, EEG_hdr, rereferenceMethod)
% REREFERENCE_EEG Re-references EEG data using specified methods
%
% This function re-references EEG data using either the 'EAR' or 'CAR' method.
%
% Inputs:
%   EEG_record - The EEG data (channels x time points)
%   EEG_hdr - Header information containing channel labels
%   rereferenceMethod - Method for re-referencing ('EAR' for linked ears or 'CAR' for common average reference)
%
% Outputs:
%   EEG_reref - The re-referenced EEG data

% Re-reference the EEG data based on the specified method
if strcmp(rereferenceMethod, 'EAR')
    % EAR (Linked Ears) Re-referencing:
    % Find the indices of the ear or mastoid channels using their labels
    channel = 1:length(EEG_hdr.label);
    A1 = channel(strcmp(EEG_hdr.label, 'A1') | strcmp(EEG_hdr.label, 'M1'));
    A2 = channel(strcmp(EEG_hdr.label, 'A2') | strcmp(EEG_hdr.label, 'M2'));
    
    % Calculate the linked ear reference by averaging the signals from A1 and A2
    linkedEar = (EEG_record(A1, :) + EEG_record(A2, :)) / 2;
    
    % Subtract the linked ear reference from each channel to re-reference the signal
    EEG_reref = EEG_record - repmat(linkedEar, size(EEG_record, 1), 1);
    
elseif strcmp(rereferenceMethod, 'CAR')
    % CAR (Common Average Reference) Re-referencing:
    % Calculate the average signal across the main 19 electrodes
    CAR = mean(EEG_record(1:19, :), 1);
    
    % Subtract the common average reference from each channel to re-reference the signal
    EEG_reref = EEG_record - repmat(CAR, size(EEG_record, 1), 1);
    
else
    % Handle invalid re-referencing method
    error('Invalid re-reference method. Choose either "EAR" or "CAR".');
end

end