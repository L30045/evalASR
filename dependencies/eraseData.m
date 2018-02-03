function EEG_erased = eraseData(EEG)
% Erase data, icaact, times, event, urevent of EEG structure to save space

EEG_erased = EEG;
EEG_erased.data = [];
EEG_erased.icaact = [];
EEG_erased.times = [];
EEG_erased.event = [];
EEG_erased.urevent = [];

end