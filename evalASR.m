%% Open EEGLAB
eeglab_path = which('eeglab.m');
if isempty(eeglab_path)
    current_path = pwd;
    addpath(['dependencies' filesep 'eeglab_current']); eeglab
    cd(current_path);
    addpath(genpath('./'));
end