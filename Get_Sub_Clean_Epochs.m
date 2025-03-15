function newStartIdx = Get_Sub_Clean_Epochs(startingIndCleanEpoch, fs, originalEpochLength, newEpochLength)
% Get_Sub_clean_Epochs calculates the starting indices of sub-epochs within larger clean epochs.
%
% Inputs:
% - startingIndCleanEpoch: A vector containing the starting indices of clean epochs in samples.
% - fs: The sampling frequency in Hz.
% - originalEpochLength: The length of the original epochs in seconds.(e.g. 30 seconds)
% - newEpochLength: The desired length of the sub-epochs in seconds. (e.g. 1 second)
%
% Output:
% - newStartIdx: A column vector containing the starting indices of each sub-epoch in samples.

    % Calculate the stopping indices for each original epoch
    % This is the last sample index of each originalEpochLength-second epoch
    stopingIndCleanEpoch = startingIndCleanEpoch + originalEpochLength * fs - 1;

    % Calculate the number of newEpochLength-second (1 second)  epochs within each original epoch
    epochRatio = originalEpochLength / newEpochLength;

    % Create a matrix of starting indices for newEpochLength-second epochs
    % Each row corresponds to an originalEpochLength-second(30 second) epoch, and each column to a newEpochLength-second (1 second) epoch
    % The bsxfun function applies element-wise subtraction to generate starting indices
    newStartIdxMatrix = bsxfun(@minus, stopingIndCleanEpoch, (fs * newEpochLength*(1:epochRatio) - 1));

    % Flatten the matrix into a single column vector
    % This converts the matrix of starting indices into a vector, listing all sub-epoch start indices
    newStartIdx = newStartIdxMatrix(:);

    % Sort the indices to ensure they are in ascending order
    % This is important if the original order was altered during matrix manipulation
    newStartIdx = sort(newStartIdx);

    % % % Alternative method with a for loop
% % % for i = 1:epochRatio
% % %     newStartIdxMatrix(i,:) = stopingIndCleanEpoch-fs*i+1; %starting from the end epoch and going down to the start of the new epoching system
% % % 
% % % end
% % % 
% % %             % Flatten the matrix into a single column vector
% % % newStartIdx = newStartIdxMatrix(:);
% % % 
% % % % Sort the indices to get them in ascending order
% % % newStartIdx = sort(newStartIdx);


end
