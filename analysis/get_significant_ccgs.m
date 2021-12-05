function [ccg_data, sig_idx] = get_significant_ccgs(ccg_data, flag)

noise_distribution2 = [ccg_data.ccg_control(:, 1:50), ccg_data.ccg_control(:, 151:201)];
noise_distribution2(2:2:end,:) = noise_distribution2(1:2:end-1,:);  
ccg_data.noise_std2 = std(noise_distribution2,0,2, 'omitnan');
ccg_data.noise_mean2 = nanmean(noise_distribution2,2);

sig_idx = (ccg_data.noise_std2>flag.sig_min_std ) & ...
        (ccg_data.peaks>(flag.sig_num_stds*ccg_data.noise_std2 + ccg_data.noise_mean2)) & ...
        (abs(ccg_data.peak_lag) <= flag.sig_max_lag);

fields = fieldnames(ccg_data);
for i = 1:length(fields)
    if ~strcmp(fields{i},'cluster') && ~strcmp(fields{i},'config') && ~strcmp(fields{i},'ccg_control')
        ccg_data.(fields{i}) = ccg_data.(fields{i})(sig_idx,:);
    end
end

