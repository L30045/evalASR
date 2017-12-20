%% Open EEGLAB
eeglab_path = which('eeglab.m');
if isempty(eeglab_path)
    current_path = pwd;
    addpath(['dependencies' filesep 'eeglab_current']); eeglab
    cd(current_path);
    addpath(genpath('./'));
end

%% Add file path
sessionPath = '/data/projects/Chunshu/NCTUStandardLevel2/session/';

%% load raw EEG data
channelToRemove = {'vehicle position'};
bandPassFilt = [0.5 100];
resamplingRate = 250;
sessNum = 1:80;
thresASR = [1, 2.5, 3, 5, 10, 20, 30, 40, 50, 70, 100, 200, 500, 1000];

for i = 1:length(sessNum)
    rawdataPath = dir([sessionPath,sprintf('%d/',sessNum(i)),'*.set']);
    rawdataPath = [rawdataPath.folder filesep rawdataPath.name];
    EEG_raw = pop_loadset(rawdataPath);
    chanLabels = {EEG_raw.chanlocs.labels};
    
    % Remove vehicle channel and reconstructed channel
    if ~isempty(EEG_raw.etc.noiseDetection.reference.interpolatedChannels)
        reconstChannel = chanLabels(EEG_raw.etc.noiseDetection.reference.interpolatedChannels);
        channelToRemove = [{'vehicle position'},reconstChannel];
    end
    EEG_raw = pop_select(EEG_raw,'nochannel',channelToRemove);
    
    % Band pass and resample
    EEG = pop_eegfiltnew(EEG, bandPassFilt(1), bandPassFilt(2), [], 0, [], 0);
    EEG = pop_resample( EEG, resamplingRate);
    
    
    
    
    
end
