function ccg_output = analyze_ccgs(data, flag)
% INPUT:
%   data: 
%   flag:
% convert spike train to 'has at least 1 spike in bin'
data.spiketrain = logical(data.spiketrain); 

% get times of interest
data_subset = data.spiketrain(:, flag.start_time:flag.end_time, :);
ccg_output = struct;

% jitter spiketrains
for j = 1:size(data_subset, 3)
    data_subset_jitter{j} = jitter(data_subset(:,:,j), flag.jit_window);
    data_subset_real{j} = data_subset(:,:,j);
end

% change to for if do not want to use multiple cpus
parfor j = 1:size(data_subset,3)
    disp("processing neuron: " + j);
    cnt = 0;
    for k = j+1:size(data_subset,3)
        % pre & post comps
        for jj = 1:2
            cnt = cnt + 1;
            if jj == 1 % first neuron is pre
                [ccgs{j}.ccg_norm(cnt,:), ccgs{j}.ccg_unnorm(cnt,:)] = xcorr_gm(data_subset_real{j}, data_subset_real{k}, flag.max_lag, flag.min_lag);
                [ccgs{j}.ccg_norm_jitter(cnt,:), ccgs{j}.ccg_unnorm_jitter(cnt,:)]  = xcorr_gm(data_subset_jitter{j}, data_subset_jitter{k}, flag.max_lag, flag.min_lag);
                ccgs{j}.neuron_id_pairs(cnt,:) = [j,k];
            else % first neuron is post
                [ccgs{j}.ccg_norm(cnt,:), ccgs{j}.ccg_unnorm(cnt,:)] = xcorr_gm(data_subset_real{k},data_subset_real{j}, flag.max_lag, flag.min_lag);
                [ccgs{j}.ccg_norm_jitter(cnt,:), ccgs{j}.ccg_unnorm_jitter(cnt,:)] = xcorr_gm(data_subset_jitter{k},data_subset_jitter{j}, flag.max_lag, flag.min_lag);
                ccgs{j}.neuron_id_pairs(cnt,:) = [k,j];
            end
        end
    end
end

% aggregate ccgs into ccg output struct
fields = fieldnames(ccgs{1});
for j = 1:length(fields)
    ccg_output.(fields{j}) = ccgs{1}.(fields{j});
end

for i = 2:size(data_subset,3)-1
    for j = 1:length(fields)
        ccg_output.(fields{j}) = [ccg_output.(fields{j}); ccgs{i}.(fields{j})];
    end
end

% ccg control is base - jitter
ccg_output.ccg_control = ccg_output.ccg_norm-ccg_output.ccg_norm_jitter;
ccg_output.ccg_control_unnorm = ccg_output.ccg_unnorm-ccg_output.ccg_unnorm_jitter;

% store input in ccg_output struct
ccg_output.flag = flag;
ccg_output.data = data;
end