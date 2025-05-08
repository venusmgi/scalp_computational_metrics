function [filterKernel] = Pick_Filter(filterType,fs)
% PICK_FILTER Selects and designs an FIR filter based on specified frequency bands.
% Rachel J. Smith (Lopouratory 2019)
% Edited by Venus (5.6.2025)
%
% This function designs FIR filters using the least-squares method for various EEG frequency
% bands ('beta', 'alpha', 'theta', 'broadband', 'delta') based on the given sampling frequency.
%
% Inputs:
%   filterType - String indicating the desired frequency band filter ('beta', 'alpha', 'theta', 'broadband', 'delta').
%   fs - Sampling frequency of the signal in Hz.
%
% Outputs:
%   filterKernel - FIR filter coefficients for the specified band.

% Compute the Nyquist frequency based on the sampling frequency
Nyq = fs / 2;

% Choose the appropriate filter design based on the specified filter type
switch lower(filterType)
    case 'beta'
        % Validate the sampling frequency for beta filter design
        if fs > 5000
            error('Sampling frequency exceeds 5000 Hz. Design a custom filter for higher frequencies.');
        end
        % Beta band filter parameters
        lowestFreq = 13; % Lowest frequency in Hz
        secPerCycle = 1/lowestFreq;
        %filter order should be at least 2 times(secPerCycle*sampleFreq)
        n = round(13*(secPerCycle*fs)); % Calculate FIR filter order
        filterBounds =  [13; 30]; % filter bounds in Hz
        % transision width should be 10-25% of the upper and lower filter bounds
        transWidth = [0.10; 0.2]; %having different transition width for different frequencies; not related to the filterBounds order
        f = [0 (1-transWidth(2))*filterBounds(1) filterBounds(1) filterBounds(2) (1+transWidth(1))*filterBounds(2) Nyq]./Nyq;
        a = [0 0 1 1 0 0]; % Amplitude response
        filterKernel = firls(n,f,a); % Design filter using least-squares method


        % Alpha filter
    case 'alpha'
        % Validate the sampling frequency for alpha filter design
        if fs > 1000
            error('Sampling frequency exceeds 1000 Hz. Design a custom filter for higher frequencies.');
        end

        % Alpha band filter parameters
        lowestFreq = 8; % Lowest frequency in Hz
        secPerCycle = 1/lowestFreq;
        %2 cycles at the lowest frequency analyzed
        n = round(16*(secPerCycle*fs)); % Calculate FIR filter order
        filterBounds =  [8; 13]; % filter bounds in Hz
        % transision width should be 10-25% of the upper and lower filter bounds
        transWidth = 0.1;
        f = [0 (1-transWidth)*filterBounds(1) filterBounds(1) filterBounds(2) (1+transWidth)*filterBounds(2) Nyq]./Nyq;
        a = [0 0 1 1 0 0]; % Amplitude response
        filterKernel = firls(n,f,a); % Design filter using least-squares method



        % Theta filter
    case 'theta'
        % Validate the sampling frequency for theta filter design
        if fs > 500
            error('Sampling frequency exceeds 500 Hz. Design a custom filter for higher frequencies.');
        end

        % Theta band filter parameters
        lowestFreq = 4;  % Lowest frequency in Hz
        secPerCycle = 1/lowestFreq;
        %2 cycles at the lowest frequency analyzed
        n = round(8*(secPerCycle*fs)); % Calculate FIR filter order; %cannot go higher than 8 because it will introduce a spike in filter response
        filterBounds = [4;8];
        transWidth = [0.15; 0.25]; %having different transition width for different frequencies; not related to the filterBounds order
        f = [0 (1-transWidth(2))*filterBounds(1) filterBounds(1) filterBounds(2) (1+transWidth(1))*filterBounds(2) Nyq]./Nyq;
        a = [0 0 1 1 0 0]; % Amplitude response
        filterKernel = firls(n,f,a); % Design filter using least-squares method


        % Broadband filter
    case 'broadband'
        % Validate the sampling frequency for broadband filter design
        if fs > 5000
            error('Sampling frequency exceeds 5000 Hz. Design a custom filter for higher frequencies.');
        end
        % Broadband filter parameters
        lowestFreq = 2; % Initial cutoff frequency in Hz
        secPerCycle = 1/lowestFreq;
        n = round(3.5*(secPerCycle*fs));% Calculate FIR filter order
        f = [0 0.25 1 55 58 Nyq]./Nyq;  % Frequency vector
        a = [0 0 1 1 0 0]; %amplitude response
        filterKernel = firls(n,f,a);



    case 'delta'
        % Validate the sampling frequency for delta filter design
        if fs > 500
            error('Sampling frequency exceeds 500 Hz. Design a custom filter for higher frequencies.');
        end
        % Delta band filter parameters

        lowestFreq = 1; % Lowest frequency in Hz
        secPerCycle = 1/lowestFreq;
        %2 cycles at the lowest frequency analyzed
        %n must be at least 2 to 5 times ( 1/lowestFreq ) * fs ;
        n = round(2*(secPerCycle*fs));
        filterBounds = [1;4];
        transWidth = 0.1; %having different transition width, to see which one will work better
        f = [0 (1-transWidth)*filterBounds(1) filterBounds(1) filterBounds(2) (1+transWidth)*filterBounds(2) Nyq]./Nyq;
        a = [0 0 1 1 0 0];
        filterKernel = fir2(n,f,a);


    otherwise
        error('Invalid filter type. Choose from: delta, theta, alpha, beta, broadband.');
end
end

