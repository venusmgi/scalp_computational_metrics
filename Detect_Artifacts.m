function [artsBool, artsStruct] = Detect_Artifacts(eeg_record,fs,stdAbove,buffer,minNumArtChans,numChans)
% Automated artifact detector based on automatic extreme value detection algorithm
% Outputs a structure with times containing artifacts
% Based on the methods of Durka et al. 2003; Moretti et al. 2003
% Adapted in Smith et al. 2021 
%
% INPUTS:
% - eeg_record:         MxN matrix containing EEG data (channels x time)    
% - fs:                 Sampling rate of the EEG
% - stdAbove:           Number of standard deviations above the EEG to be
%                       considered an artifact (default: 7.5)
% - buffer:             Time in seconds before and after a detected artifact
%                       that will be marked as artifact (default: 0.9s)
% - minNumArtChans:     Minimum number of channels with excessive artifact
%                       to count and index as artifactual (default: 1)
% - numChans:           Number of channels to include in artifact
%                       detection; channels from index 1 up to this number
%                       will be included, e.g. numchans = 19 includes
%                       channels 1:19
% 
% OUTPUTS:
% - artsStruct:         Structure containing the times and type of
%                       artifacts seen (artifacts.times and artifacts.general)
% - artsBool:           boolean array marking the indices of the EEG signal
%                       that are artifactual; 1 means artifact and 0 means clean
% 
%
% Previously called: "artifact_detector_RJS_041119.m"
% Made by Rachel J. Smith (Lopouratory, 2019)
%
% V.1.0: Removed 19 channel restriction, treats all rows as EEG
% V.1.1: Uploaded onto Lopouratory Github, 2022)
% Updated by Venus Mostaghimi (Lopouratory, 2025)

% Minimum threshold; artifact must be > 200 uV regardless of threshold
% parameter
minThresh = 200;

% Take only the first numChans channels 
eeg_record = eeg_record(1:numChans,:);

% Filter the EEG so the appearnace matches the clinical viewer
% Create filters - these were compared to Nihon Kohden viewer in Aug 2015
tau = 0.1;  % time constant in Nihon Kohden
cutoff = 1/(2*pi*tau); % cutoff frequency is a function of time constant
[b1,a1] = butter(1, cutoff/(fs/2), 'high');
[b2,a2] = butter(3, 40/(fs/2), 'low');

% Note: To view frequency response of filter you can use freqz as below
% freqz(b1,a1,[],fs)

% Filter the data - Nihon Kohden does not use zero-phase filters
% Note: In spring 2025, we compared this to the equivalent zero-phase
% filters and found that the time shift when using the causal filters was
% negligible. Also, when using the causal filters, the algorithm was
% slightly more sensitive to artifacts (the zero-phase filters dampened the
% signals more), so we kept the causal filters despite being unconventional.
eeg_filt = filter(b1,a1,eeg_record',[],1);
eeg_filt = -filter(b2,a2,eeg_filt,[],1); % POSITIVE DOWN
eeg_filt = eeg_filt';

% Calc mean and standard deviation over time to determine threshold
meanEEG = mean(eeg_filt,2);
stdEEG = std(eeg_filt,0,2);

% Ensure that artifacts are bigger than the minimum threshold, after
% removing the mean; if values in stdEEG are too small, make them bigger
stdEEG(stdEEG<(minThresh/stdAbove)) = minThresh/stdAbove;

possArts = abs(eeg_filt - meanEEG) > stdAbove*stdEEG;
% possArts = abs(eeg_filt - repmat(meanEEG,1,size(eeg_filt,2))) > repmat(stdAbove,1,size(eeg_filt,2)).*stdEEG;


% Get logical vector for artifacts occurring in minimum # of channels
sumArtsVec = sum(possArts,1);
artsVec = (sumArtsVec > minNumArtChans-1);


% Get logical vector for impedance check artifacts
diff_EEGrec = diff(eeg_record,1,2);
impChecks = diff_EEGrec == 0;
sumImpCheck = sum(impChecks,1);
impCheckVec = sumImpCheck>8;

% Combine general artifacts with impedance check artifacts
impCheckVec = [false impCheckVec]; % pad to match the size, no artifact at the begining
artsBool = artsVec|impCheckVec;


% Pad around each artifact to ensure you get the entire event (buffer = # of seconds)
indicesArts = find(artsBool);
for k=1:size(indicesArts,2)
    if indicesArts(k)-buffer*fs<=0 % artifact too close to the beginning
        artsBool(1:indicesArts(k)+buffer*fs) = true;
    elseif indicesArts(k)+buffer*fs>=size(eeg_record,2) % artifact too close to the end
        artsBool(indicesArts(k)-buffer*fs:end) = true;
    else  
        artsBool(indicesArts(k)-buffer*fs:indicesArts(k)+buffer*fs) = true;
    end
end


    
% If the user requests the artsStruct output, which turns the artifactual
% indices into start and end times and includes info on the artifact type
if nargout > 1
    % Convert logical to times that can be read into the artifact_marking code
    fullIndArts = find(artsBool); % all artifact indices
    if ~isempty(fullIndArts)
        % Determine number of artifacts
        diffVec = diff(fullIndArts); % difference will be 1 during an artifact, ~=1 in between successive artifacts
        largeDiff = find(diffVec~=1);  % find boundaries between artifacts
        numArtifacts = length(largeDiff) + 1; % # of artifacts is number of boundaries + 1
        
        % Create artTimes matrix; loop through all artifacts and calculate
        % start and end times in seconds
        artTimes = zeros(numArtifacts,2);
        artTimes(1,1) = fullIndArts(1)./fs; % beginning of the first artifact
        for m = 1:size(largeDiff,2)
            artTimes(m,2) = fullIndArts(largeDiff(m))./fs; % end of the current artifact
            artTimes(m+1,1) = fullIndArts(largeDiff(m)+1)./fs; % beginning of the next artifact
        end
        artTimes(size(largeDiff,2)+1,2) = fullIndArts(end)./fs; % end of the last artifact
    else
        % If no artifacts, return zeros
        artTimes = zeros(0,2);
    end
    
    % Create output structure
    artsStruct.times = artTimes;  % artifact start/end times in seconds
    artsStruct.impedance = zeros(1,size(artTimes,1)); % artifacts due to impedance checks
    
    % Mark any artifact that came from impedence check
    if ~isempty(artTimes)
        for i = 1:size(artTimes,1)
            startSample = round(artTimes(i,1)*fs); % start index of artifact i
            endSample = round(artTimes(i,2)*fs); % end index of artifact i
            if any(impCheckVec(startSample:endSample))
                artsStruct.impedance(i) = 1;
            end
        end
    end
    artsStruct.general = double(~artsStruct.impedance);

end