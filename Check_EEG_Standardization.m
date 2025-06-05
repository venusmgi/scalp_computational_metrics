function Check_EEG_Standardization(desiredChannelOrder, headerEEG)
% CHECK_EEG_STANDARDIZATION Verifies the EEG channel order and reference names.
%
% This function checks whether the data contained in `headerEEG` follows the
% specified channel order defined by `desiredChannelOrder`. It ensures that the EEG
% channels are in the correct sequence and that reference channels (A1/M1 and A2/M2)
% are appropriately labeled. If the channel order is not standardized, it raises an
% error with suggestions on how to address the issue.
%
% Inputs:
%   desiredChannelOrder - A cell array specifying the desired standard order of EEG channels.
%   headerEEG - A structure containing EEG header information, such as channel labels.
%
% Outputs:
%   None. An error is raised if the EEG channel order is not standardized.

% Determine the number of channels specified by the desired order
    numChannels = length(desiredChannelOrder);
    
% Check if the ear/mastoid channel names match the desired naming convention
actualEarChannels = ismember({'A1', 'A2', 'M1', 'M2'}, headerEEG.label);
expectedEarChannels = ismember({'A1', 'A2', 'M1', 'M2'}, desiredChannelOrder);

if ~isequal(actualEarChannels, expectedEarChannels)
    warning(['The ear/mastoid channel names differ from the desired standard.', ...
        ' Check the channel name to ensure accuracy in channel labeling. Ignoring this difference for now.']);
    % Handle potential mismatches between expected and actual reference names.
    % For example, if labels include A1/A2 but desired order includes
    % M1/M2, change desired order to match labels
    if all(ismember({'A1', 'A2'}, headerEEG.label)) && all(ismember({'M1', 'M2'}, desiredChannelOrder))
        desiredChannelOrder(ismember(desiredChannelOrder, {'M1', 'M2'})) = {'A1', 'A2'};
    elseif all(ismember({'M1', 'M2'}, headerEEG.label)) && all(ismember({'A1', 'A2'}, desiredChannelOrder))
        desiredChannelOrder(ismember(desiredChannelOrder, {'A1', 'A2'})) = {'M1', 'M2'};
    end
end

% Compare the current channel order against the desired order
isChannelOrderCorrect = strcmp(headerEEG.label(1:numChannels), desiredChannelOrder);

%% if we don't want to have the same order and we just need to make sure that a spesific channel exists
% Compare the current channel order against the desired order
% for i = 1:length(headerEEG.label)
%     isChannelOrderCorrect(i) = ismember(headerEEG.label{i}, desiredChannelOrder);
% end
%%


% Verify whether all channels are in the desired order
if ~all(isChannelOrderCorrect)
    error(['The EEG channel order is not standardized.', ...
        ' Please use the "Standardize_EEG_Channel_Order" function to correct the channel order and labels.']);
end
