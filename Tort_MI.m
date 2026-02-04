function [modulationIndex, correctnessFlag] = Tort_MI(amplitudeValues, phaseValues, numPhaseBins)
% Tort_tMI - Calculate Modulation Index using Tort's method for Phase-Amplitude Coupling
%
% DESCRIPTION:
%   Computes the Modulation Index (MI) as described by Tort et al. (2010).
%   This metric quantifies phase-amplitude coupling by measuring how much the
%   amplitude distribution across phase bins deviates from a uniform distribution.
%   The MI is based on the Kullback-Leibler (KL) divergence between the observed
%   amplitude distribution and a uniform distribution, normalized by the maximum
%   possible entropy.
%
% SYNTAX:
%   [modulationIndex, correctnessFlag] = Tort_tMI(amplitudeValues, phaseValues, numPhaseBins)
%
% INPUTS:
%   amplitudeValues - Vector of instantaneous amplitude values from high-frequency
%                     signal (e.g., gamma band amplitude from Hilbert transform)
%   phaseValues     - Vector of instantaneous phase values from low-frequency
%                     signal in degrees (range: -180 to 180)
%   numPhaseBins    - Number of phase bins for discretization (typically 18,
%                     corresponding to 20-degree bins as suggested in literature)
%
% OUTPUTS:
%   modulationIndex  - Modulation Index (range: 0 to 1)
%                      0 = no coupling (uniform amplitude distribution)
%                      1 = maximum coupling (all amplitude in one phase bin)
%   correctnessFlag  - Quality check flag
%                      0 = calculation correct (distribution sums to 1)
%                      1 = potential numerical issue detected
%
% ALGORITHM:
%   1. Bin phases into equal-width bins
%   2. Calculate mean amplitude for each phase bin
%   3. Normalize to create probability distribution P(φ)
%   4. Calculate KL divergence: D_KL = log(N) - H(P)
%      where H(P) is Shannon entropy and N is number of bins
%   5. Normalize: MI = D_KL / log(N)
%
% EXAMPLE:
%   % Calculate MI for gamma amplitude coupled to delta phase
%   hilbertGamma = hilbert(gammaSignal);
%   amplitude = abs(hilbertGamma);
%   hilbertDelta = hilbert(deltaSignal);
%   phase = angle(hilbertDelta) * 180/pi;  % Convert to degrees
%   numBins = 18;
%   [MI, isCorrect] = Tort_tMI(amplitude, phase, numBins);
%
% REFERENCE:
%   Tort, A. B., Komorowski, R., Eichenbaum, H., & Kopell, N. (2010).
%   Measuring phase-amplitude coupling between neuronal oscillations of
%   different frequencies. Journal of Neurophysiology, 104(2), 1195-1210.
%
% SEE ALSO:
%   canoltyMVL, getPAC
%
% AUTHOR: Venus Mostaghimi

%% Phase binning
% Create phase bins from -180 to 180 degrees
% Literature suggests using 18 bins of 20 degrees each
phaseBinEdges = linspace(-180, 180, numPhaseBins + 1);

% Assign each phase value to its corresponding bin
binnedPhaseIndices = discretize(phaseValues, phaseBinEdges);

%% Assign amplitudes to phase bins
% Create matrix where each row represents a phase bin
amplitudeBins = nan(numPhaseBins, length(phaseValues));

for binIdx = 1:numPhaseBins
    % Store amplitude values that fall within current phase bin
    amplitudeBins(binIdx, binnedPhaseIndices == binIdx) = ...
        amplitudeValues(binnedPhaseIndices == binIdx);
end

%% Calculate mean amplitude for each phase bin
% Compute mean amplitude per bin, ignoring NaN values
meanAmplitudePerBin = mean(amplitudeBins, 2, 'omitnan');

% Replace NaN with zero for bins with no data points
meanAmplitudePerBin(isnan(meanAmplitudePerBin)) = 0;

%% Normalize amplitude distribution to create probability distribution
% Sum of mean amplitudes across all bins
totalMeanAmplitude = sum(meanAmplitudePerBin);

% Create normalized probability distribution P(φ)
normalizedDistribution = meanAmplitudePerBin / totalMeanAmplitude;

%% Verify distribution validity
% Check if probability distribution sums to 1 (within floating-point precision)
distributionSum = sum(normalizedDistribution);
distributionSum = single(distributionSum);

correctnessFlag = 0;
if distributionSum ~= 1
    correctnessFlag = 1;  % Flag potential numerical issue
end

%% Calculate Kullback-Leibler (KL) divergence
epsilon = 1e-6;  % Small value to avoid log(0) numerical instability

% Calculate Shannon entropy: H(P) = -sum(P(i) * log(P(i)))
entropyTerms = zeros(numPhaseBins, 1);
for binIdx = 1:numPhaseBins
    % Add epsilon to prevent log(0)
    entropyTerms(binIdx) = normalizedDistribution(binIdx) * ...
                           log(normalizedDistribution(binIdx) + epsilon);
end

% Total Shannon entropy
shannonEntropy = -sum(entropyTerms);

% Calculate KL divergence from uniform distribution
% D_KL = log(N) - H(P), where N is number of bins
uniformEntropy = log(numPhaseBins);  % Maximum possible entropy
klDivergence = uniformEntropy - shannonEntropy;

%% Calculate normalized Modulation Index
% Normalize by maximum possible KL divergence to get value between 0 and 1
modulationIndex = klDivergence / uniformEntropy;

end