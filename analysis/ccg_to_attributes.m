function ccg_data = ccg_to_attributes(ccg_output, flag)
ccg_types = ["ccg"];

for ccg_type = 1:length(ccg_types)
    ccg_data_t = struct;
    ccg_data_t.ccg_control = ccg_output.(strcat(ccg_types(ccg_type),'_control'));
    ccg_data_t.ccg_jitter = ccg_output.(strcat(ccg_types(ccg_type),'_norm_jitter'));
    ccg_data_t.ccg_norm = ccg_output.(strcat(ccg_types(ccg_type),'_norm'));
    ccg_data_t.ccg_unnorm = ccg_output.(strcat(ccg_types(ccg_type),'_unnorm'));

    ccgs = ccg_data_t.ccg_control(:,91:111);
    pair_ids = ccg_output.neuron_id_pairs;
    data_included = ccg_output.data;
    
    ccg_data_t.session = ccg_output.session;
    ccg_data_t.cluster.Cluster_session = data_included.Cluster_session;
    
    maxlag = flag.max_pt_lag;
    
    [peaks, peak_lag] = max(ccgs,[],2);
    [troughs, trough_lag] = min(ccgs, [], 2);
    [peaks, lpeak_lag] = max(ccg_data_t.ccg_control,[],2);
    [troughs, ltrough_lag] = min(ccg_data_t.ccg_control,[],2);
    ccg_data_t.config.cont_vars = ["peaks", "troughs", "peak_lag","trough_lag", "area", "peak_trough_diff", "peak_trough_lag_diff","pair_distance","peak_width", "trough_width"];
    ccg_data_t.config.cat_vars = ["pre_id", "post_id", "pre_ct", "post_ct", "pre_cl", "post_cl", "pre_depth", "post_depth", "clust", "syn_peak", "syn_trough"];
    ccg_data_t.config.flag = flag;    
    ccg_data_t.peaks = peaks;
    ccg_data_t.troughs = troughs;
    
    
    ccg_data_t.peak_lag = 100-lpeak_lag+1;
    ccg_data_t.trough_lag = 100-ltrough_lag+1;
        
    ccg_data_t.peak_trough_diff =  ccg_data_t.peaks- ccg_data_t.troughs;
    ccg_data_t.peak_trough_lag_diff = ccg_data_t.peak_lag-ccg_data_t.trough_lag;   
    ccg_data_t.syn_peak = ccg_data_t.peak_lag == 0; 
    ccg_data_t.syn_trough = ccg_data_t.trough_lag == 0; 
    
    %% computing r-orientation, r-eye, and r-sf for each pair of neurons
    pre_id = pair_ids(:,1);
    post_id = pair_ids(:,2);
    
    pair_positions = data_included.Cluster_celldepth(pair_ids);
    ccg_data_t.pair_distance  = abs(pair_positions(:,1)-pair_positions(:,2));
    
    ccg_data_t.pre_id = pre_id;
    ccg_data_t.post_id = post_id;
    
    if isfield(data_included, 'layer_centers') && flag.debug_mode == 0
        ccg_data_t.config.layer_centers = data_included.layer_centers;
        ccg_data_t.pre_depth_scale = data_included.Cluster_celldepth_scale(ccg_data_t.pre_id);
        ccg_data_t.post_depth_scale =  data_included.Cluster_celldepth_scale(ccg_data_t.post_id);
    end
    
    ccg_data_t.pre_ct = data_included.Cluster_celltype(ccg_data_t.pre_id);
    ccg_data_t.post_ct = data_included.Cluster_celltype(ccg_data_t.post_id);
    ccg_data_t.config.ct_labels = ["AS", "FS", "RM", "RL"];

    if flag.split_layers
        ccg_data_t.config.cl_labels = flip([ "WM", "6","5", "4cB", "4cA", "4a/b", "2/3"]);
    else
        ccg_data_t.config.cl_labels = flip(["WM", "5/6", "4c", "4a/b", "2/3"]);
    end
    
    ccg_data_t.config.sc_labels = ["complex", "simple"];
    ccg_data_t.cluster.Cluster_celllayer = data_included.Cluster_celllayer;
    ccg_data_t.cluster.Cluster_celldepth = data_included.Cluster_celldepth;
    ccg_data_t.cluster.Cluster_celltype = data_included.Cluster_celltype;
    ccg_data_t.cluster.Cluster_MI_max = data_included.Cluster_MI_max;
    ccg_data_t.cluster.Cluster_simpcomp = data_included.Cluster_simpcomp;
    
    % time to first spike
    cluster_mean_time = [];
    for i = 1:length(data_included.spiketrain)
        st = logical(data_included.spiketrain{i});
        for j = 1:size(st,3)
            st_neu = st(:,:,j);
            [ms, is] = max((st_neu(:,101:350)), [], 2);
            is(ms==0) = NaN;
            cluster_mean_time = [cluster_mean_time; nanmedian(is)];
        end
    end
    ccg_data_t.cluster.time_to_fs = cluster_mean_time;
    
    
    ccg_data_t.pre_sc = data_included.Cluster_simpcomp(ccg_data_t.pre_id);
    ccg_data_t.post_sc = data_included.Cluster_simpcomp(ccg_data_t.post_id);
    
    ccg_data_t.pre_mi = data_included.Cluster_MI_max(ccg_data_t.pre_id);
    ccg_data_t.post_mi = data_included.Cluster_MI_max(ccg_data_t.post_id);
    
    ccg_data_t.pre_cl = data_included.Cluster_celllayer(ccg_data_t.pre_id);
    ccg_data_t.post_cl = data_included.Cluster_celllayer(ccg_data_t.post_id);
    
    ccg_data_t.pre_depth = data_included.Cluster_celldepth(ccg_data_t.pre_id);
    ccg_data_t.post_depth = data_included.Cluster_celldepth(ccg_data_t.post_id);
    
    eyes = 0:2;
    for j = 1:length(eyes)
        cnt = 1;
        for i = 1:length(ccg_output.data.spiketrain)
            eye_idx = (floor((ccg_output.data.conditionID{i}-1)/144) == eyes(j))';
            st = ccg_output.data.spiketrain{i};
            mean_eye_resp(cnt:cnt+size(st,3)-1,j) = mean(st(eye_idx, 100:end, :), [1,2]);
            cnt = cnt + size(st,3);
        end
    end
    
    orientations = 0:35;
    for j = 1:length(orientations)
        cnt = 1;
        for i = 1:length(ccg_output.data.spiketrain)
            orientation_idx = mod(ccg_output.data.conditionID{i}-1,36) == orientations(j);
            st = ccg_output.data.spiketrain{i};
            mean_ori_resp(cnt:cnt+size(st,3)-1,j) = mean(st(orientation_idx, 100:end, :), [1,2]);
            cnt = cnt + size(st,3);
        end
    end
    
    freqs = 0:3;
    for j = 1:length(freqs)
        cnt = 1;
        for i = 1:length(ccg_output.data.spiketrain)
            freq_idx = (floor(mod(ccg_output.data.conditionID{i}-1,144)/36) == freqs(j))';
            st = ccg_output.data.spiketrain{i};
            mean_freq_resp(cnt:cnt+size(st,3)-1,j) = mean(st(freq_idx, 100:end, :), [1,2]);
            cnt = cnt + size(st,3);
        end
    end
        
    for i = 1:length(ccg_data_t.pre_id)
        ccg_data_t.r_ori(i) = corr(mean_ori_resp(pre_id(i),:)', mean_ori_resp(post_id(i),:)');
        ccg_data_t.r_freq(i) = corr(mean_freq_resp(pre_id(i),:)', mean_freq_resp(post_id(i),:)');
        ccg_data_t.r_eye(i) = corr(mean_eye_resp(pre_id(i),:)', mean_eye_resp(post_id(i),:)');
    end
    
    ccg_data_t.cluster.mean_ori_tune = mean_ori_resp;
    ccg_data_t.cluster.mean_freq_tune = mean_freq_resp;
    ccg_data_t.cluster.mean_eye_tune = mean_eye_resp;

    %% ccg_exclusions
    noise_distribution2 = [ccg_data_t.ccg_control(:, 1:50), ccg_data_t.ccg_control(:, 151:201)];
    ccg_data_t.area = nansum(ccgs,2);
    ccg_data_t.noise_std2 = std(noise_distribution2,[],2, 'omitnan');
    ccg_data_t.noise_mean2 = nanmean(noise_distribution2,2);
    
    %% ccg peak and trough width
    ccg_data_t.peak_width = zeros(1,size(ccgs,1));
    ccg_data_t.trough_width = zeros(1,size(ccgs, 1));
    
    for i = 1:size(ccgs, 1)
        ccg = ccgs(i,:);
        % peak width
        ccg_left = ccg(1:peak_lag(i));
        
        ccg_right = ccg(peak_lag(i):end);
        left_ind = find((ccg_left<peaks(i)-abs(peaks(i)/2)), 1, 'last')+1;
        right_ind = find((ccg_right<peaks(i)-abs(peaks(i)/2)), 1, 'first')+(peak_lag(i)-1)-1;
        if isempty(left_ind)
            left_ind = 1;
        end
        if isempty(right_ind)
            right_ind = maxlag+1;
        end
        ccg_data_t.peak_width(i) = right_ind-left_ind;
        
        % trough width
        ccg_left = ccg(1:trough_lag(i));
        left_ind = find((ccg_left>troughs(i)+abs(troughs(i)/2)), 1, 'last')+1;
                
        ccg_right = ccg(trough_lag(i):end);
        right_ind = find((ccg_right>troughs(i) + abs(troughs(i)/2)), 1, 'first')-1+(trough_lag(i)-1);
        
        if isempty(left_ind)
            left_ind = 1;
        end
        if isempty(right_ind)
            right_ind = maxlag+1;
        end
        ccg_data_t.trough_width(i) = right_ind-left_ind;
    end
    
    if flag.large_output
    else
        ccg_data_t.ccg_control = [];
    end
    ccg_data_t.ccgs = ccgs;
    
    % iterate over output struct and make everything verticle
    fields = fieldnames(ccg_data_t);
    for i = 1:length(fields)
        if ~strcmp(fields{i},'cluster') && ~strcmp(fields{i},'config')
            if size(ccg_data_t.(fields{i}), 1) == 1
                ccg_data_t.(fields{i}) = ccg_data_t.(fields{i})';
            end
            %ccg_data_t.(fields{i}) = single(ccg_data_t.(fields{i}));
        end
    end
    ccg_data.(strcat(ccg_types(ccg_type))) = ccg_data_t;
end
end