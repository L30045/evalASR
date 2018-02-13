%% Open EEGLAB
restoredefaultpath
eeglab_path = which('eeglab.m');
if isempty(eeglab_path)
    current_path = pwd;
    addpath(['dependencies' filesep 'eeglab_current']); eeglab
    cd(current_path);
    addpath(genpath('./'));
end

%% Add file path
sessionPath = '/data/projects/yuan/lane-keeping task dataset with chanlocs/';
saveSessPath = '/data/projects/yuan/evalASR/simplified dataset/';
saveParasPath = '/home/yuan/Documents/evalASR/parameters/';

%% load raw EEG data
channelToRemove = {'vehicle position'};
bandPassFilt = [0.5 100];
resamplingRate = 250;
thresASR = [1, 2.5, 3, 5, 10, 20, 30, 40, 50, 70, 100, 200, 500, 1000];

%% begin main pipeline
for i = 1:1
    % Structure to save all the parameters
    if ~exist([saveParasPath,sprintf('paras_s%d.mat',i)],'file')
        paras = initParas();
        paras.sessNum = i;
    else
        load([saveParasPath,sprintf('paras_s%d.mat',i)]);
    end
    
    rawdataPath = dir([sessionPath,sprintf('s%d/',i),'*.set']);
    rawdataPath = [rawdataPath.folder filesep rawdataPath.name];
    EEG_raw = pop_loadset(rawdataPath);
    
    % Remove vehicle channel and reconstructed channel
    channelToRemove = {'vehicle position'};
    
    
    EEG_raw = pop_select(EEG_raw,'nochannel',channelToRemove);
    nCh = EEG_raw.nbchan;
    rankData = rank(double(EEG_raw.data(:,1:1e5)));
    paras.etc.rankData = rankData;
    paras.etc.channelRemoved = channelToRemove;
    paras.etc.chanlocs = EEG_raw.chanlocs;
    
    % Band pass and resample
    EEG_raw = pop_eegfiltnew(EEG_raw, bandPassFilt(1), bandPassFilt(2), [], 0, [], 0);
    EEG_raw = pop_resample( EEG_raw, resamplingRate);
    
    %%
    if rankData == nCh
    
    %% Raw section
    % compute var for each channels
    paras.channel_variance.var_raw_all = var(EEG_raw.data,0,2);
    
    % compute ICs for raw data
%     if rankData ~= nCh
%         [para.icaweights.W_raw, para.icasphere.S_raw] = runica(double(EEG_raw.data),'lrate',0.001,'extended',1,'pca',rankData);
%     else
%         [para.icaweights.W_raw, para.icasphere.S_raw] = runica(double(EEG_raw.data),'lrate',0.001,'extended',1);
%     end
    [paras.icaweights.W_raw, paras.icasphere.S_raw] = runica(double(EEG_raw.data),'extended',1);
    
    % compute marginal entropy
    paras.entropy.ent_raw = getent2(EEG_raw.data);
    icaact = paras.icaweights.W_raw * paras.icasphere.S_raw * EEG_raw.data;
    paras.entropy.ent_rawICA = getent2(icaact);
    
    %% Clean section
    [EEG_clean, clean_mask, paras.etc.clean_info] = getClean(EEG_raw,i);
    % compute var for raw data
    paras.channel_variance.var_raw_cSection = var(EEG_raw.data(:,clean_mask),0,2);
    paras.channel_variance.var_raw_bSection = var(EEG_raw.data(:,~clean_mask),0,2);
    paras.channel_variance.var_raw_all = var(EEG_raw.data,0,2);
    
    % compute ICs for clean section 
%     if rankData ~= nCh
%         [para.icaweights.W_clean, para.icasphere.S_clean] = runica(double(EEG_clean.data),'lrate', 0.001,'extended',1,'pca',rankData);
%     else
%         [para.icaweights.W_clean, para.icasphere.S_clean] = runica(double(EEG_clean.data),'lrate', 0.001,'extended',1);
%     end
    [paras.icaweights.W_clean, paras.icasphere.S_clean] = runica(double(EEG_clean.data),'extended',1);    
    
    % compute marginal entropy
    paras.entropy.ent_clean = getent2(EEG_clean.data);
    icaact = paras.icaweights.W_clean * paras.icasphere.S_clean * EEG_clean.data;
    paras.entropy.ent_cleanICA = getent2(icaact);
    
    % compute clean ICs activites on raw data
    paras.icaact.icaact_cleanIC_cSection = var(paras.icaweights.W_clean * paras.icasphere.S_clean * EEG_clean.data,0,2);
    paras.icaact.icaact_cleanIC_bSection = var(paras.icaweights.W_clean * paras.icasphere.S_clean * EEG_raw.data(:,~clean_mask).data,0,2);
    paras.icaact.icaact_cleanIC_all = var(paras.icaweights.W_clean * paras.icasphere.S_clean * EEG_raw.data,0,2);
    
    % compute raw ICs activites on raw data
    paras.icaact.icaact_rawIC_cSection = var(paras.icaweights.W_raw * paras.icasphere.S_raw * EEG_raw.data(:,clean_mask),0,2);
    paras.icaact.icaact_rawIC_bSection = var(paras.icaweights.W_raw * paras.icasphere.S_raw * EEG_raw.data(:,~clean_mask),0,2);
    paras.icaact.icaact_rawIC_all = var(paras.icaweights.W_raw * paras.icasphere.S_raw * EEG_raw.data,0,2);
    
    %% save simplified EEG dataset
    savepath = sprintf([saveSessPath,'s%d/'],i);
    mkdir(savepath);
    EEG_raw = eraseData(EEG_raw);
    EEG_clean = eraseData(EEG_clean);
    eval(sprintf('pop_saveset(EEG_raw, ''filename'', ''s%d_resampled.set'',''filepath'',savepath);', i));
    eval(sprintf('pop_saveset(EEG_clean, ''filename'', ''s%d_cleanSection.set'',''filepath'',savepath);', i));
    
    %% ASR section
    % disable high pass filter, bad channel removal and window rejection for ASR
    thresFlatChannel = -1;
    highPassBand = -1;
    thresPoorCorrChannel = -1;
    thresLineNoiseChannel = -1;
    thresWindow = -1;
    
    % spatial_filter for raw ICs
    spatial_filter_raw = paras.icaweights.W_raw * paras.icasphere.S_raw;
    % spatial_filter for clean ICs
    spatial_filter_clean = paras.icaweights.W_clean * paras.icasphere.S_clean;
    
    % para for parfor
    var_asr_clean = zeros(nCh,length(thresASR));
    var_asr_bad = zeros(nCh,length(thresASR));
    var_asr_all = zeros(nCh,length(thresASR));
    act_rawIC_clean = zeros(nCh,length(thresASR));
    act_rawIC_bad = zeros(nCh,length(thresASR));
    act_rawIC_all = zeros(nCh,length(thresASR));
    act_cleanIC_clean = zeros(nCh,length(thresASR));
    act_cleanIC_bad = zeros(nCh,length(thresASR));
    act_cleanIC_all = zeros(nCh,length(thresASR));
    act_asr1IC_clean = zeros(nCh,length(thresASR));
    act_asr1IC_bad = zeros(nCh,length(thresASR));
    act_asr1IC_all = zeros(nCh,length(thresASR));
    w_asr = zeros(nCh, nCh, length(thresASR));
    s_asr = zeros(nCh, nCh, length(thresASR));
    ent_asr = zeros(nCh, length(thresASR));
    ent_asrICA = zeros(nCh, length(thresASR));
    
    %% compute thresASR = 1
    % ASR with different threshold
    EEG_asr = clean_rawdata(EEG_raw, thresFlatChannel, highPassBand, thresPoorCorrChannel, thresLineNoiseChannel, thresASR(1), thresWindow);

    % run ICA
%     if rankData ~= nCh
%         [w_1, s_1] = runica(double(EEG_asr.data),'lrate', 0.001,'extended',1,'pca',rankData);
%     else
%         [w_1, s_1] = runica(double(EEG_asr.data),'lrate', 0.001,'extended',1);
%     end
    [w_1, s_1] = runica(double(EEG_asr.data),'extended',1);
    w_asr(:,:,1) = w_1;
    s_asr(:,:,1) = s_1;
    
    % spatial_filter for thresASR = 1 ICs
    spatial_filter_asr1 = w_1 * s_1;
    
    % compute var for each channels
    var_asr_all(:,1) = var(EEG_asr.data,0,2);
    var_asr_clean(:,1) = var(EEG_asr.data(:,clean_mask),0,2);
    var_asr_bad(:,1) = var(EEG_asr.data(:,~clean_mask),0,2);

    % compute raw ICs activities
    act_rawIC_clean(:,1) = var(spatial_filter_raw*EEG_asr.data(:,clean_mask),0,2);
    act_rawIC_bad(:,1) = var(spatial_filter_raw*EEG_asr.data(:,~clean_mask),0,2);
    act_rawIC_all(:,1) = var(spatial_filter_raw*EEG_asr.data,0,2);

    % compute clean ICs activites
    act_cleanIC_clean(:,1) = var(spatial_filter_clean*EEG_asr.data(:,clean_mask),0,2);
    act_cleanIC_bad(:,1) = var(spatial_filter_clean*EEG_asr.data(:,~clean_mask),0,2);
    act_cleanIC_all(:,1) = var(spatial_filter_clean*EEG_asr.data,0,2);

    % compute thresASR = 1 ICs activites
    act_asr1IC_clean(:,1) = var(spatial_filter_asr1*EEG_asr.data(:,clean_mask),0,2);
    act_asr1IC_bad(:,1) = var(spatial_filter_asr1*EEG_asr.data(:,~clean_mask),0,2);
    act_asr1IC_all(:,1) = var(spatial_filter_asr1*EEG_asr.data,0,2);
    
    % compute marginal entropy
    ent_asr(:,1) = getent2(EEG_asr.data);
    ent_asrICA(:,1) = getent2(spatial_filter_asr1*EEG_asr.data);
    
    % save simplified EEG dataset
    EEG_asr = eraseData(EEG_asr);
    filename = sprintf('s%d_ASR%.1f.set',i, thresASR(1));
    pop_saveset(EEG_asr,'filename', filename,'filepath',savepath);
    
    %% compute thresASR    
    % begin parpool
%     parpool(length(thresASR));

    for j = 2:length(thresASR)
        % ASR with different threshold
        EEG_asr = clean_rawdata(EEG_raw, thresFlatChannel, highPassBand, thresPoorCorrChannel, thresLineNoiseChannel, thresASR(j), thresWindow);
        
        % run ICA
        [w_temp, s_temp] = runica(double(EEG_asr.data),'extended',1);
        w_asr(:,:,j) = w_temp;
        s_asr(:,:,j) = s_temp;
        
        % compute var for each channels
        var_asr_all(:,j) = var(EEG_asr.data,0,2);
        var_asr_clean(:,j) = var(EEG_asr.data(:,clean_mask),0,2);
        var_asr_bad(:,j) = var(EEG_asr.data(:,~clean_mask),0,2);
        
        % compute raw ICs activities
        act_rawIC_clean(:,j) = var(spatial_filter_raw*EEG_asr.data(:,clean_mask),0,2);
        act_rawIC_bad(:,j) = var(spatial_filter_raw*EEG_asr.data(:,~clean_mask),0,2);
        act_rawIC_all(:,j) = var(spatial_filter_raw*EEG_asr.data,0,2);
        
        % compute clean ICs activites
        act_cleanIC_clean(:,j) = var(spatial_filter_clean*EEG_asr.data(:,clean_mask),0,2);
        act_cleanIC_bad(:,j) = var(spatial_filter_clean*EEG_asr.data(:,~clean_mask),0,2);
        act_cleanIC_all(:,j) = var(spatial_filter_clean*EEG_asr.data,0,2);
        
        % compute thresASR = 1 ICs activites
        act_asr1IC_clean(:,j) = var(spatial_filter_asr1*EEG_asr.data(:,clean_mask),0,2);
        act_asr1IC_bad(:,j) = var(spatial_filter_asr1*EEG_asr.data(:,~clean_mask),0,2);
        act_asr1IC_all(:,j) = var(spatial_filter_asr1*EEG_asr.data,0,2);
     
        % compute marginal entropy
        ent_asr(:,j) = getent2(EEG_asr.data);
        ent_asrICA(:,j) = getent2(w_temp*s_temp*EEG_asr.data);
        
        % save simplified EEG dataset
        EEG_asr = eraseData(EEG_asr);
        filename = sprintf('s%d_ASR%.1f.set',i, thresASR(j));
        pop_saveset(EEG_asr,'filename', filename,'filepath',savepath);
    end
    
    % close parpool
%     delete(gcp);
    %% save para to para structure    
    for k = 1:length(thresASR)
        eval(sprintf('para.ASR.var.var_cSection_%.1f = var_asr_clean(:,k)',thresASR(k)));
        eval(sprintf('para.ASR.var.var_bSection_%.1f = var_asr_bad(:,k)',thresASR(k)));
        eval(sprintf('para.ASR.var.var_all_%.1f = var_asr_all(:,k)',thresASR(k)));
        eval(sprintf('para.ASR.cleanIC.icaact_cleanIC_cSection_%.1f = act_cleanIC_clean(:,k)',thresASR(k)));
        eval(sprintf('para.ASR.cleanIC.icaact_cleanIC_bSection_%.1f = act_cleanIC_bad(:,k)',thresASR(k)));
        eval(sprintf('para.ASR.cleanIC.icaact_cleanIC_all_%.1f = act_cleanIC_all(:,k)',thresASR(k)));
        eval(sprintf('para.ASR.rawIC.icaact_rawIC_cSection_%.1f = act_rawIC_clean(:,k)',thresASR(k)));
        eval(sprintf('para.ASR.rawIC.icaact_rawIC_bSection_%.1f = act_rawIC_bad(:,k)',thresASR(k)));
        eval(sprintf('para.ASR.rawIC.icaact_rawIC_all_%.1f = act_rawIC_all(:,k)',thresASR(k)));
        eval(sprintf('para.ASR.icaweights.W_%.1f = w_asr(:,:,k)',thresASR(k)));
        eval(sprintf('para.ASR.icasphere.S_%.1f = s_asr(:,:,k)',thresASR(k)));
        eval(sprintf('para.entropy.ent_asr_%.1f = ent_asr(:,k)',thresASR(k)));
        eval(sprintf('para.entropy.ent_asrICA_%.1f = ent_asrICA(:,k)',thresASR(k)));
    end
    
    else
        notTestyet = [notTestyet, i];
    end
    
    %% save parameters
    eval(sprintf('para_s%02d = para;',i));
    if ~exist('parameters/para_info.mat', 'file')
        save('parameters/para_info.mat',sprintf('para_s%02d',i));
    else
        save('parameters/para_info.mat',sprintf('para_s%02d',i),'-append');
    end
    
end
