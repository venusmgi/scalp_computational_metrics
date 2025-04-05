% Check_EEG_Standardization
% Verifies that the EEG channel order matches the desired standard order.
%
% This function checks whether the EEG data in `headerEEG` follows the specified
% channel order in `desiredChannelOrder`. It ensures that the channels are in the
% correct sequence and that reference channels (A1/M1 and A2/M2) are correctly labeled.
% If the channel order is not standardized, an error is raised with a suggestion
% to use a specific function for standardization.
%
% Inputs:
%   desiredChannelOrder - A cell array specifying the desired order of EEG channels.
%   headerEEG - A structure containing EEG header information, including channel labels.
%
% Outputs:
%   None. The function raises an error if the channel order is not standardized.

function Check_EEG_Standardization(desiredChannelOrder, headerEEG)

    % Initialize the array to track correct channel order
    lenChanOrder = length(desiredChannelOrder);
    correctChannOrder = nan(lenChanOrder, 1);

    % Check if the first 19 channels match the desired order
    for i = 1:19
        correctChannOrder(i, 1) = isequal(headerEEG.label(i), desiredChannelOrder(i));
    end

    % Check if the 20th channel is either "A1" or "M1"
    correctChannOrder(20, 1) = isequal(headerEEG.label(20), "A1") || isequal(headerEEG.label(20), "M1");

    % Check if the 21st channel is either "A2" or "M2"
    correctChannOrder(21, 1) = isequal(headerEEG.label(21), "A2") || isequal(headerEEG.label(21), "M2");

    % Verify that all channels are in the correct order
    if sum(correctChannOrder) ~= lenChanOrder
        error("The EEG channel order is not standardized. Use the function 'Get_desired_channel_order_output_excel.m' to standardize the EDF and convert it to a MAT file, then upload the EEG.");
    end
end