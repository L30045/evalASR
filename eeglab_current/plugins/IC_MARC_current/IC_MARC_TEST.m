% test IC_MARC
addpath('/home/yuan/Documents/eeglab_current');
eeglab
EEG = pop_loadset('Eeglab_data.set');
EEG = eeg_checkset(EEG);
EEG = pop_chanedit(EEG,'load',{'/home/yuan/Documents/eeglab_current/sample_data/eeglab_chan32.locs' 'filetype' 'autodetect'});
EEG = eeg_checkset(EEG);
EEG = pop_eegfiltnew(EEG,[],1,424,true,[],0);
EEG = eeg_checkset(EEG);
EEG = pop_selectevent(EEG,'position',2,'deleteevents','on');
EEG.setname= 'Square, Postion 2';
EEG = eeg_checkset(EEG);
EEG = pop_runica(EEG,'extended',1,'interupt','on');
% pop_saveset(EEG,'/home/yuan/Downloads/Eeglab_data.set');
pop_selectcomps(EEG,1:32);
EEG = pop_epoch(EEG,{'square'},[-1 2], 'newname',...
    'Square, Postion 2 epochs', 'epochinfo','yes');
EEG = pop_rmbase(EEG,[-1000 0]);
eeglab redraw