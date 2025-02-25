function [autoArts] = get_automatedArtifacts_EEG(eeg_record,fs,stdAbove,buffer,channelsInvolved)
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
% - channelsInvolved:   Number of channels with excessive artifact to be
%                       counted as artifact (default =1)
% 
% OUTPUTS:
% - autoArts:           Structure containing the times and type of
%                       artifacts seen (artifacts.times and artifacts.general)
%
% Previously called: "artifact_detector_RJS_041119.m"
% Made by Rachel J. Smith (Lopouratory, 2019)
%
% V.1.0: Removed 19 channel restriction, treats all rows as EEG
% V.1.1: Uploaded onto Lopouratory Github, 2022)

% Removed 19 channel restriction
% % Take only the first 19 channels 
% eeg_record = eeg_record(1:19,:);

    % filter data the same way they do to clinically view
    % Create filters - these were compared to Nihon Kohden viewer in Aug 2015
    tau = 0.1;  % time constant in Nihon Kohden
    cutoff = 1/(2*pi*tau); % cutoff frequency is a function of time constant
    [b,a] = butter(1, cutoff/(fs/2), 'high');
    [b2,a2] = butter(3, 40/(fs/2), 'low');

    % Filter the data - Nihon Kohden does not use zero-phase filters
    eeg_filt = filter(b,a,eeg_record',[],1);
    eeg_filt = -filter(b2,a2,eeg_filt,[],1); % POSITIVE DOWN
    eeg_filt = eeg_filt';

    meanEEG = mean(eeg_filt,2);
    stdEEG = std(eeg_filt,0,2);

    % set an absolute threshold that an artifact must be over 200 uV. 
    stdEEG(stdEEG<(200/stdAbove)) = 200/stdAbove;

    possArts = zeros(size(eeg_filt));
    for j=1:size(eeg_filt,1)
        possArts(j,:) = abs(eeg_filt(j,:)-meanEEG(j))>stdAbove*stdEEG(j);
    end

    % get logical vector for any artifacts
    sumArtsVec = sum(possArts,1);
    artsVec = sumArtsVec>channelsInvolved-1;

    % pad around the artifact to ensure you get the entire event (buffer = # of seconds)
    indicesArts = find(artsVec);
    for k=1:size(indicesArts,2)
        if indicesArts(k)-buffer*fs<=0
            artsVec(1:indicesArts(k)+buffer*fs) = 1;
        elseif indicesArts(k)+buffer*fs>=size(eeg_record,2)
            artsVec(indicesArts(k)-buffer*fs:end) = 1;
        else
            artsVec(indicesArts(k)-buffer*fs:indicesArts(k)+buffer*fs) = 1;
        end
    end

    % get it into times that can be read into the artifact_marking code
    % difference between indices vector will show which have a difference more
    % than 1
    fullIndArts = find(artsVec);
    diffVec = diff(fullIndArts);
    largeDiff = find(diffVec~=1); %find when the artifact sample is non-consecutive
    artTimes = [];
    if ~isempty(fullIndArts)
        artTimes = zeros(1,2);
        artTimes(1,1) = fullIndArts(1)./fs; 
        for m = 1:size(largeDiff,2)
            artTimes(m,2) = fullIndArts(largeDiff(m))./fs;
            artTimes(m+1,1) = fullIndArts(largeDiff(m)+1)./fs;
        end
        artTimes(size(largeDiff,2)+1,2) = fullIndArts(end)./fs;
    end

    % add artifacts from impedance checks
    diff_EEGrec = diff(eeg_record,1,2);
    diffIsZero = diff_EEGrec == 0;
    impChecks = double(diffIsZero);
    sumImpCheck = sum(impChecks);
    impCheckVec = sumImpCheck>8;

    impCheckArts = find(impCheckVec==1);

    newArts = [];
    if ~isempty(impCheckArts)
        diffVec2 = diff(impCheckArts);
        largeDiff2 = find(diffVec2~=1); %find when the artifact sample is non-consecutive

        newArts = zeros(1,2);
        newArts(1,1) = impCheckArts(1)./fs; 
        for p = 1:size(largeDiff2,2)
         newArts(p,2) = impCheckArts(largeDiff2(p))./fs;
         newArts(p+1,1) = impCheckArts(largeDiff2(p)+1)./fs;
        end
        newArts(size(largeDiff2,2)+1,2) = impCheckArts(end)./fs;
    end
    artTimes = [artTimes; newArts];

    autoArts.times = artTimes;
    autoArts.general = ones(1,size(artTimes,1));
    autoArts.eye = zeros(1,size(artTimes,1));
    autoArts.ear = zeros(1,size(artTimes,1));

end