function [EEG_reref] = Rereference_EEG(EEG_record, EEG_hdr, rereferenceMethod)
% REREFERENCE_EEG Re-references EEG data using specified methods.
%
% This function re-references EEG data using one of several methods:
% 'EAR' for linked ears, 'CAR' for common average reference, 'bipolar' for bipolar referencing,
% or a specific channel name for single-channel referencing.
%
% Inputs:
%   EEG_record - Matrix of EEG data (channels x time points).
%   EEG_hdr - Header information containing channel labels.
%   rereferenceMethod - Method for re-referencing:
%       'EAR' - Linked ears re-referencing.
%       'CAR' - Common average reference.
%       'bipolar' - Bipolar re-referencing.
%       Channel name (e.g., 'Cz') for single-channel referencing.
%
% Outputs:
%   EEG_reref - The re-referenced EEG data.


% Define the list of channel names for potential single-channel referencing
channelList =  {'Fp1','Fp2','F3','F4','C3','C4','P3','P4','O1',...
    'O2','F7','F8','T3','T4','T5','T6','Fz','Cz','Pz','M1','M2','Eye1','Eye2','EKG1','EKG2','EMG1','EMG2'};

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
    
elseif strcmp(rereferenceMethod, 'bipolar')
    % Bipolar Re-referencing:
    % Calculate differences between pairs of channels for bipolar montage
    % bipolar_pairs = {'Fp1-F3', 'F3-C3', 'C3-P3', 'P3-O1',...
    %              'Fp1-F7','F7-T3', 'T3-T5','T5-O1',...
    %              'Fp2-F4', 'F4-C4','C4-P4', 'P4-O2', ...
    %              'Fp2-F8','F8-T4', 'T4-T6', 'T6-O2'};

    Fp1 = Get_Channel_Loc ('Fp1',EEG_hdr);
    Fp2 = Get_Channel_Loc ('Fp2',EEG_hdr);
    F3 = Get_Channel_Loc ('F3',EEG_hdr);
    F4 =Get_Channel_Loc ('F4',EEG_hdr);
    C3 = Get_Channel_Loc ('C3',EEG_hdr);
    C4 = Get_Channel_Loc ('C4',EEG_hdr);
    P3 = Get_Channel_Loc ('P3',EEG_hdr);
    P4 = Get_Channel_Loc ('P4',EEG_hdr);
    O1 = Get_Channel_Loc ('O1',EEG_hdr);
    O2 = Get_Channel_Loc ('O2',EEG_hdr);
    F7 = Get_Channel_Loc ('F7',EEG_hdr);
    F8 =Get_Channel_Loc ('F8',EEG_hdr);
    T3 = Get_Channel_Loc ('T3',EEG_hdr);
    T4 = Get_Channel_Loc ('T4',EEG_hdr);
    T5 = Get_Channel_Loc ('T5',EEG_hdr);
    T6 = Get_Channel_Loc ('T6',EEG_hdr);



    EEG_reref = nan(16,length(EEG_record));
    EEG_reref(1,:) = EEG_record(Fp1,:)-EEG_record(F3,:);
    EEG_reref(2,:) = EEG_record(F3,:)-EEG_record(C3,:);
    EEG_reref(3,:) = EEG_record(C3,:)-EEG_record(P3,:);
    EEG_reref(4,:) = EEG_record(P3,:)-EEG_record(O1,:);
    
    EEG_reref(5,:) = EEG_record(Fp1,:)-EEG_record(F7,:);
    EEG_reref(6,:) = EEG_record(F7,:)-EEG_record(T3,:);
    EEG_reref(7,:) = EEG_record(T3,:)-EEG_record(T5,:);
    EEG_reref(8,:) = EEG_record(T5,:)-EEG_record(O1,:);
    
    
    EEG_reref(9,:) = EEG_record(Fp2,:)-EEG_record(F4,:);
    EEG_reref(10,:) = EEG_record(F4,:)-EEG_record(C4,:);
    EEG_reref(11,:) = EEG_record(C4,:)-EEG_record(P4,:);
    EEG_reref(12,:) = EEG_record(P4,:)-EEG_record(O2,:);
    
    EEG_reref(13,:) = EEG_record(Fp2,:)-EEG_record(F8,:);
    EEG_reref(14,:) = EEG_record(F8,:)-EEG_record(T4,:);
    EEG_reref(15,:) = EEG_record(T4,:)-EEG_record(T6,:);
    EEG_reref(16,:) = EEG_record(T6,:)-EEG_record(O2,:);


elseif contains(rereferenceMethod, channelList) 
    % Single-Channel Re-referencing:
    % Check if onl a single channel is chosen
    if sum(contains(rereferenceMethod, channelList)) >1
        error ('Choose only one channel for re-referencing.')
    end
    % Find the reference channel index
    referenceChannel = strcmp(rereferenceMethod,channelList);
    % Subtract the reference channel from all channels
    EEG_reref = EEG_record - EEG_record(referenceChannel,:);

else
    % Handle invalid re-referencing method
    error('Invalid re-reference method. Choose either "EAR", "CAR", "bipolar", or a channel name such as "Cz".');
end

end

%% Helper Functions

% GET_CHANNEL_LOC Finds the index of a specified channel label in the EEG header
function channLoc = Get_Channel_Loc (channel_label,EEG_hdr)
    channLoc = strcmp(EEG_hdr.label, channel_label);
end