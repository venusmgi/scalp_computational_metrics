function [autoArts,artsVec] = get_automatedArtifacts_EEG(eeg_record,fs,stdAbove,buffer,minNumArtChannels)
% Automated artifact detector based on automatic extreme value detection algorithm
% Outputs a structure with times containing artifacts
% Based on the methods of Durka et al. 2003; Moretti et al. 2003
% Adapted in Smith et al. 2021 
%
% INPUTS:
% - eeg_record:         MxN matrix containing EEG data (channels x time).    
% - fs:                 Sampling rate of the EEG
% - stdAbove:           Number of standard deviations above the EEG to be considered an artifact (default: 7.5)
% - buffer:             Time in seconds before and after detected an artifact to
%                       remove (default 0.9s)
% - channelsInvolved:   Minimum number of channels with excessive artifact to be
%                       counted as artifact (default =1)
% 
% OUTPUTS:
% - autoArts:           Structure containing the times and type of
%                       artifacts seen (artifacts.times and artifacts.general)
% - artsVec:            logical array, marking the indices of the EEG signal that are artifactual,
%                       1 means artifact and 0 means clean index
% 
%
% Previously called: "artifact_detector_RJS_041119.m"
% Made by Rachel J. Smith (Lopouratory, 2019)
%
% V.1.0: Removed 19 channel restriction, treats all rows as EEG
% V.1.1: Uploaded onto Lopouratory Github, 2022)

% Removed 19 channel restriction
% Take only the first 19 channels 
eeg_record = eeg_record(1:19,:);

% filter data the same way they do to clinically view
% Create filters - these were compared to Nihon Kohden viewer in Aug 2015
tau = 0.1;  % time constant in Nihon Kohden
cutoff = 1/(2*pi*tau); % cutoff frequency is a function of time constant
[b1,a1] = butter(1, cutoff/(fs/2), 'high');

% freqz(b,a,[],fs)
[b2,a2] = butter(3, 40/(fs/2), 'low');
% freqz(b2,a2,[],fs)

% Filter the data - Nihon Kohden does not use zero-phase filters
eeg_filt = filter(b1,a1,eeg_record',[],1);
eeg_filt = -filter(b2,a2,eeg_filt,[],1); % POSITIVE DOWN
eeg_filt = eeg_filt';

meanEEG = mean(eeg_filt,2);
stdEEG = std(eeg_filt,0,2);


% set an absolute threshold that an artifact must be over 200 uV. 
stdEEG(stdEEG<(200/stdAbove)) = 200/stdAbove;

possArts = false(size(eeg_filt));
for j=1:size(eeg_filt,1)
    possArts(j,:) = abs(eeg_filt(j,:)-meanEEG(j))>stdAbove*stdEEG(j);
end

% get logical vector for any artifacts
sumArtsVec = sum(possArts,1);
artsVec = sumArtsVec>minNumArtChannels-1;

% get logical vector for impedance check artifacts
diff_EEGrec = diff(eeg_record,1,2);
impChecks = diff_EEGrec == 0;
sumImpCheck = sum(impChecks);
impCheckVec = sumImpCheck>8;

impCheckVec = [false;impCheckVec]; % pad to match the size, nor artifact at the begining
totalArts = or(artsVec,impCheckVec); % combine general artifacts and impedance artifacts

% pad around the artifact to ensure you get the entire event (buffer = # of seconds)
indicesArts = find(totalArts);
if ~isempty(indicesArts)
    for k = 1:length(indicesArts)
        startIdx = max(1, indicesArts(k)-buffer*fs);
        endIdx = min(length(indicesArts),indicesArts(k)+buffer*fs);
        artsVec(startIdx:endIdx) = true;

    end
end


% get it into times that can be read into the artifact_marking code
% difference between indices vector will show which have a difference more
% than 1
fullIndArts = find(artsVec);
if ~isempty(fullIndArts)
    % Pre-allocate artTimes
    diffVec = diff(fullIndArts);
    largeDiff = find(diffVec~=1);  %find when the artifact sample is non-consecutive
    numArtifacts = length(largeDiff) + 1;
    artTimes = zeros(numArtifacts,2);
    artTimes(1,1) = fullIndArts(1)./fs; %begining of the first artifact
    for m = 1:size(largeDiff,2)
        artTimes(m,2) = fullIndArts(largeDiff(m))./fs; %ending of the current artifact
        artTimes(m+1,1) = fullIndArts(largeDiff(m)+1)./fs; %begining of the next artifact
    end
    artTimes(size(largeDiff,2)+1,2) = fullIndArts(end)./fs;%ending of the last artifact
else
    artTimes = zeros(0,2);
end

% Create output structure
autoArts.times = artTimes;
autoArts.general = ones(1,size(artTimes,1));
% autoArts.eye = zeros(1,size(artTimes,1));
% autoArts.ear = zeros(1,size(artTimes,1));
autoArts.impedance = zeros(1,size(artTimes,1));

%marking any artifact that came from impedence check
if ~isempty(artTimes)
    for i = 1:length(artTimes)
        startSmaple = round(artTimes(i,1)*fs);
        endSample = round(artTimes(i,2)*fs);
        if any(impCheckVec(startSmaple:endSample))
            autoArts.impedance(i) = 1;
        end
    end
end

end