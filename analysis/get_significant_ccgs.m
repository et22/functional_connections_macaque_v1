function [ccg_data, sig_idx] = get_significant_ccgs(ccg_data, flag)

sig_idx = (ccg_data.noise_std2>flag.sig_min_std ) & ...
        (ccg_data.peaks>(flag.sig_num_stds*ccg_data.noise_std2 + ccg_data.noise_mean2)) & ...
        (abs(ccg_data.peak_lag) <= flag.sig_max_lag);

fields = fieldnames(ccg_data);
for i = 1:length(fields)
    if ~strcmp(fields{i},'cluster') && ~strcmp(fields{i},'config') && ~strcmp(fields{i},'ccg_control')
        ccg_data.(fields{i}) = ccg_data.(fields{i})(sig_idx,:);
    end
end

