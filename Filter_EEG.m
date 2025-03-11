function [filteredEEG] = Filter_EEG(EEG_record, fs, filterType)
% FILTER_EEG Filters EEG data using a specified filter type
%
% This function applies a specified filter to EEG data.
%
% Inputs:
%   EEG_record - The EEG data (channels x time points)
%   fs - Sampling frequency of the EEG data
%   filterType - Type of filter to apply (e.g., 'lowpass', 'highpass', etc.)
%
% Outputs:
%   filteredEEG - The filtered EEG data (channels x time points)

% Filter the EEG signal
% Pick_Filter is assumed to be a custom function that returns the filter coefficients
filter = Pick_Filter(filterType, fs);

% Apply the filter using zero-phase filtering
% filtfilt applies the filter in both forward and reverse directions to eliminate phase distortion
% EEG_record' transposes the data so that filtfilt operates along the columns (time points)
% The result is transposed back to maintain the original orientation (channels x time points)
filteredEEG = filtfilt(filter, 1, EEG_record')';

end