%% Open EEGLAB
restoredefaultpath
eeglab_path = which('eeglab.m');
if isempty(eeglab_path)
    addpath('eeglab_current/'); eeglab
end
addpath(genpath('dependencies/'));
%% Add file path
sessionPath = '/data/projects/yuan/lane-keeping task dataset with chanlocs/';
saveSessPath = '/data/projects/yuan/evalASR/simplified dataset/';
saveParasPath = '/home/yuan/Documents/evalASR/parameters/';

%% process parameters setting
channelToRemove = {'vehicle position'};
bandPassFilt = [0.5 100];
resamplingRate = 250;
thresASR = [1, 2.5, 3, 5, 10, 20, 30, 40, 50, 70, 100, 200, 500, 1000];
parFlag = 0;

%% begin main pipeline
for i = 1:28
    % Structure to save all the parameters
    if ~exist([saveParasPath,sprintf('paras_s%02d.mat',i)],'file')
        paras = initParas();
        paras.sessNum = i;
    else
        load([saveParasPath,sprintf('paras_s%02d.mat',i)]);
    end
    paras.etc.ASR_threshold = thresASR;
    
    rawdataPath = dir([sessionPath,sprintf('s%d/',i),'*.set']);
    rawdataPath = [rawdataPath.folder filesep rawdataPath.name];
    EEG_raw = pop_loadset(rawdataPath);

    % Band pass and resample
    EEG_resampled = pop_eegfiltnew(EEG_raw, bandPassFilt(1), bandPassFilt(2), [], 0, [], 0);
    EEG_resampled = pop_resample( EEG_resampled, resamplingRate);
    paras.preprocessing_info.ori_samplingRate = EEG_resampled.srate;
    paras.preprocessing_info.resamplingRate = resamplingRate;
    paras.preprocessing_info.bandPassFilt = bandPassFilt;
    
    % Remove vehicle channel and reference channel 
    channelToRemove = {'vehicle position','A1','A2'};    
    EEG_resampled = pop_select(EEG_resampled,'nochannel',channelToRemove);
    ori_chanlocs = EEG_resampled.chanlocs;
    % Add chanlocs
    if exist([rawdataPath(1:end-3),'xyz'],'file')
        EEG_resampled.chanlocs = readlocs([rawdataPath(1:end-3),'xyz'],'filetype','xyz');
    elseif exist([rawdataPath(1:end-3),'DAT'],'file')
        EEG_resampled.chanlocs = readlocs([rawdataPath(1:end-3),'DAT'],'filetype','DAT');
    else
        disp('No EEG channel location found');
    end
    paras.etc.chanlocs = EEG_resampled.chanlocs;
    
    % bad channel removal
    thresFlatChannel = 5;
    highPassBand = -1;
    thresPoorCorrChannel = 0.8;
%     thresPoorCorrChannel = 0.85;
    thresLineNoiseChannel = 4;
    thresWindow = -1;
    EEG_resampled = clean_rawdata(EEG_resampled, thresFlatChannel, highPassBand, thresPoorCorrChannel, thresLineNoiseChannel, -1, thresWindow);
    if isfield(EEG_resampled.etc,'clean_channel_mask')
        paras.preprocessing_info.removedChannel = [channelToRemove, {ori_chanlocs(~EEG_resampled.etc.clean_channel_mask).labels}];
    else
        paras.preprocessing_info.removedChannel = channelToRemove;
    end
    nCh = EEG_resampled.nbchan;
    rankData = rank(EEG_resampled.data(:,1:1e5));
    
    %% Get original channel variance
    % compute var for each channels in whole section
    paras.ori_Ch_var.ori_var_wholeSec = var(EEG_resampled.data,0,2);
    
    % get clean section
    [EEG_clean, mask, paras.etc.clean_windows_info.clean_portion] = getClean(EEG_resampled);
    paras.etc.clean_windows_info.mask = mask;
    
    % compute var for each channels in clean section
    paras.ori_Ch_var.ori_var_cleanSec = var(EEG_clean.data,0,2);
    
    %% Get ICs from raw data in whole section and clean section
%     if rankData ~= nCh
%         [para.icaweights.W_raw, para.icasphere.S_raw] = runica(double(EEG_raw.data),'lrate',0.001,'extended',1,'pca',rankData);
%     else
%         [para.icaweights.W_raw, para.icasphere.S_raw] = runica(double(EEG_raw.data),'lrate',0.001,'extended',1);
%     end

    % ICs from whole section
    EEG_resampled = pop_runica(EEG_resampled, 'icatype', 'runica', 'extended',1);
    [w_whole, s_whole] = deal(EEG_resampled.icaweights, EEG_resampled.icasphere);
    paras.ori_ICA_result.whole_Section.icaweights = w_whole;
    paras.ori_ICA_result.whole_Section.icasphere = s_whole;
    
    % ICs from clean section
    EEG_clean = pop_runica(EEG_clean, 'icatype', 'runica', 'extended',1);
    [w_clean, s_clean] = deal(EEG_clean.icaweights, EEG_clean.icasphere);
    paras.ori_ICA_result.clean_Section.icaweights = w_clean;
    paras.ori_ICA_result.clean_Section.icasphere = s_clean;
    
    % whole section ICs' activities in whole section and clean section
    paras.ori_ICA_result.whole_Section.ICact_ori = var(w_whole*s_whole*EEG_resampled.data,0,2);
    paras.ori_ICA_result.whole_Section.ICact_ori_clean = var(w_clean*s_clean*EEG_resampled.data,0,2);
    
    % clean section ICs' activities in whole section and clean section
    paras.ori_ICA_result.clean_Section.ICact_ori = var(w_whole*s_whole*EEG_clean.data,0,2);
    paras.ori_ICA_result.clean_Section.ICact_ori_clean = var(w_clean*s_clean*EEG_clean.data,0,2);
    
    %% save resampled EEG dataset
    savepath = sprintf([saveSessPath,'s%d/'],i);
    mkdir(savepath);
    eval(sprintf('pop_saveset(EEG_resampled, ''filename'', ''s%d_resampled.set'',''filepath'',savepath);', i));
    eval(sprintf('pop_saveset(EEG_clean, ''filename'', ''s%d_resampled_cleanSec.set'',''filepath'',savepath);', i));    

    %% save parameters
    eval(sprintf('paras_s%02d = paras;',i));
    save([saveParasPath, sprintf('paras_s%02d',i)],sprintf('paras_s%02d',i));
    
    
end
