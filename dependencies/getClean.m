function [clean_section, clean_mask, clean_info] = getClean(EEG, sessNum)
% Get clean mask for dataset
% Input:
%       EEG     data to clean
%       sessNum session index
% Output:
%       clean_section
%         EEG structure with concatenated clean windows
%       clean_info
%         a structure contain:
%           1. session index
%           2. clean_mask
%           3. portion of data preserved

ref_maxbadchannels = 0.075;
ref_tolerances = [-3.5 5.5];
ref_wndlen = 1;
eval(sprintf('clean_info_s%02d = struct([]);',sessNum));

clean_info = struct;
clean_info.sessNum = sessNum;

[clean_section ,clean_mask] = clean_windows(EEG,ref_maxbadchannels,ref_tolerances,ref_wndlen);
clean_info.clean_mask = clean_mask;

portion_preserved = sum(clean_mask)/length(clean_mask);
clean_info.portion_preserved = portion_preserved;

% eval(sprintf('clean_info_s%02d = clean_info;',sessNum));

% if ~exist('parameters/clean_windows_info.mat', 'file')
%     save('parameters/clean_windows_info.mat',sprintf('clean_info_s%02d',sessNum));
% else
%     save('parameters/clean_windows_info.mat',sprintf('clean_info_s%02d',sessNum),'-append');
% end

end