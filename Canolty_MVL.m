function [mvl] = Canolty_MVL(amplitudeValues, phaseValues)

% Canolty_MVL - Calculate Mean Vector Length for Phase-Amplitude Coupling
%
% DESCRIPTION:
%   Computes the Mean Vector Length (MVL) as described by Canolty et al. (2006).
%   This metric quantifies the strength of phase-amplitude coupling by measuring
%   the consistency of amplitude modulation across phase angles. The MVL represents
%   the normalized magnitude of the complex-valued mean of amplitude-weighted
%   phase vectors.
%
% SYNTAX:
%   [mvl] = Canolty_MVL(amplitudeValues, phaseValues)
%
% INPUTS:
%   amplitudeValues - Vector of instantaneous amplitude values from high-frequency
%                     signal (e.g., gamma band amplitude from Hilbert transform)
%   phaseValues     - Vector of instantaneous phase values from low-frequency
%                     signal in radians (e.g., delta phase from Hilbert transform)
%                     Note: If input is in degrees, convert to radians first
%
% OUTPUTS:
%   mvl - Mean Vector Length (normalized, range: 0 to max amplitude)
%         Higher values indicate stronger phase-amplitude coupling
%         Value of 0 indicates no coupling (uniform phase distribution)
%
% ALGORITHM:
%   MVL = |mean(A * e^(i*φ))| = |sum(A * e^(i*φ))| / N
%   where:
%     A = amplitude envelope of high-frequency signal
%     φ = phase of low-frequency signal
%     N = number of samples
%     | | denotes absolute value (magnitude of complex vector)
%
% EXAMPLE:
%   % Calculate MVL for gamma amplitude coupled to delta phase
%   hilbertGamma = hilbert(gammaSignal);
%   amplitude = abs(hilbertGamma);
%   hilbertDelta = hilbert(deltaSignal);
%   phase = angle(hilbertDelta);  % Phase in radians
%   mvl = Canolty_MVL(amplitude, phase);
%
% REFERENCE:
%   Canolty, R. T., Edwards, E., Dalal, S. S., Soltani, M., Nagarajan, S. S.,
%   Kirsch, H. E., ... & Knight, R. T. (2006). High gamma power is phase-locked
%   to theta oscillations in human neocortex. Science, 313(5793), 1626-1628.
%
% SEE ALSO:
%   tortMI, Calc_PAC
%
% AUTHOR: Venus Mostaghimi

% Calculate complex-valued composite signal: amplitude-weighted phase vectors
compositeSignal = amplitudeValues(:) .* exp(1i * phaseValues(:));



% Calculate Mean Vector Length (normalized by number of samples)
numSamples = length(amplitudeValues);
mvl = abs(sum(compositeSignal)) / numSamples;





end