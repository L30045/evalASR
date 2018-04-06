% process ASR without doing ICA
% in order to calculate IC_raw and IC_clean activity

% Process ASR with different threshold and save parameters
%% Open EEGLAB
restoredefaultpath
eeglab_path = which('eeglab.m');
if isempty(eeglab_path)
    addpath('eeglab_current/'); eeglab
end
addpath(genpath('dependencies/'));
%% Add file path
saveSessPath = '/data/projects/yuan/evalASR/simplified dataset/';
saveParasPath = '/home/yuan/Documents/evalASR/parameters/';

thresASR = [1, 2.5, 3, 5, 10, 20, 30, 40, 50, 70, 100, 200, 500, 1000];

%% ASR process
for i = 1:5
    % Structure to save all the parameters
    load([saveParasPath,sprintf('paras_s%02d.mat',i)]);
    eval(sprintf('paras = paras_s%02d;',i));
    if isfield(paras.ASR_result,'clean_Sectioin')
        paras.ASR_result.clean_Section = paras.ASR_result.whole_Section;
        paras.ASR_result = rmfield(paras.ASR_result,'clean_Sectioin');
    end
    
    EEG_resampled = pop_loadset([saveSessPath,sprintf('s%d/s%d_resampled.set',i,i)]);
    EEG_clean = pop_loadset([saveSessPath,sprintf('s%d/s%d_resampled_cleanSec.set',i,i)]);
    
    savepath = sprintf([saveSessPath,'s%d/'],i); 

    %% ASR section setting
    % disable high pass filter, bad channel removal and window rejection for ASR
    thresFlatChannel = -1;
    highPassBand = -1;
    thresPoorCorrChannel = -1;
    thresLineNoiseChannel = -1;
    thresWindow = -1;
    
    w_whole = paras.ori_ICA_result.whole_Section.icaweights;
    s_whole = paras.ori_ICA_result.whole_Section.icasphere;
    w_clean = paras.ori_ICA_result.clean_Section.icaweights;
    s_clean = paras.ori_ICA_result.clean_Section.icasphere;
    
    % spatial_filter whole ICs
    spatial_filter_whole = w_whole * s_whole;
    % spatial_filter clean ICs
    spatial_filter_clean = w_clean * s_clean;
    % clean mask
    mask = paras.etc.clean_windows_info.mask;
    
    %% ASR process
    ori_struct = paras.ASR_result;
    for j = 1:length(thresASR)
        if j==1
            previous_result = [];
        else
            previous_result = paras.ASR_result;
        end
        paras.ASR_result = ori_struct;
        
        % ASR with different threshold
        EEG_asr = clean_rawdata(EEG_resampled, thresFlatChannel, highPassBand, thresPoorCorrChannel, thresLineNoiseChannel, thresASR(j), thresWindow);
        paras.cutoff = thresASR(j);
%         % run ICA
%         EEG_asr = pop_runica(EEG_asr, 'icatype', 'runica', 'extended',1);
%         [w_asr, s_asr] = deal(EEG_asr.icaweights, EEG_asr.icasphere);
%         paras.ASR_result.icaweights = w_asr;
%         paras.ASR_result.icasphere = s_asr;
%         % spatial_filter ASR ICs
%         spatial_filter_ASR = w_asr * s_asr;

        % data modified
        paras.whole_Section.dataModified = nnz(EEG_asr.data - EEG_resampled.data)/numel(EEG_resampled.data);
        paras.clean_Section.dataModified = nnz(EEG_asr.data(:,mask) - EEG_clean.data)/numel(EEG_clean.data);

        % compute var for each channels
        paras.whole_Section.var = var(EEG_asr.data,0,2);
        paras.clean_Section.var = var(EEG_asr.data(:,mask),0,2);

        % compute raw ICs activities
        paras.whole_Section.ICact_ori = var(spatial_filter_whole*EEG_asr.data,0,2);
        paras.clean_Section.ICact_ori = var(spatial_filter_whole*EEG_asr.data(:,mask),0,2);

        % compute clean ICs activites
        paras.whole_Section.ICact_ori_clean = var(spatial_filter_clean*EEG_asr.data,0,2);
        paras.clean_Section.ICact_ori_clean = var(spatial_filter_clean*EEG_asr.data(:,mask),0,2);

%         % compute ASR ICs activites
%         paras.ASR_result.whole_Section.ICact_asr = var(spatial_filter_ASR*EEG_asr.data,0,2);
%         paras.ASR_result.clean_Section.ICact_asr = var(spatial_filter_ASR*EEG_asr.data(:,mask),0,2);

%         % save simplified EEG dataset
%         EEG_asr = eraseData(EEG_asr);
%         filename = sprintf('s%d_ASR%.1f.set',i, thresASR(j));
%         pop_saveset(EEG_asr,'filename', filename,'filepath',savepath);

       paras.ASR_result = [previous_result, paras.ASR_result]; 
    end
    
    
    %% save parameters
    eval(sprintf('paras_s%02d = paras;',i));
    save([saveParasPath, sprintf('paras_s%02d',i)],sprintf('paras_s%02d',i),'-append');
    
    
end
