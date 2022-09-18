%% This script generates plots for the supplement regarding the distribution of clusters 
% across layers. 

%% load cluster output
flag = config();

for group_idx = 1:flag.group_cnt
    clearvars -except flag group_idx
    
    load(flag.cluster_output(group_idx), 'clust_data');
    load(flag.postproc_output(group_idx), 'ccg_data');
    
    if ~exist(flag.figure_dir, 'dir')
        mkdir(flag.figure_dir)
    end
    
    flag.figure_dir = flag.figure_dir + "figS6_";

    %% subsetting ccg_data based on significance
    ccg_data = ccg_data.ccg;

    [ccg_data, sig_idx] = get_significant_ccgs(ccg_data, flag);
    ccg_data.cl_labels = ["2/3", "4a/b", "4c\alpha", "4c\beta", "5", "6", "WM"];
    ccg_data.sc_labels = ["comp.", "simp."];
    ccg_data.ct_labels = ["AS", "FS", "RM", "RL"];
    %% set up plotting
    cmap = [14,121,178;243,146,55;191,19,99;244,159,188;6,214,160]./255;%[linspace(247, 198, clust_data.num_clusters)', linspace(219,121, clust_data.num_clusters)', linspace(187,88,clust_data.num_clusters)']./255;
    cmap = [cmap; 1-cmap];
    cmap=cmap(1:clust_data.num_clusters,:);
    [empty_label{1:clust_data.num_clusters}] = deal('');
    clust_data.labels = sort_clusters_by_lag(clust_data.num_clusters, ccg_data.ccgs, clust_data.labels);
    
    rng(1);
    
    
    %% generate final_idx, random direction for CCGs to use for figures
    % generate random half
    rng(2); % seed rng for replicability
    inc_id = rand(1,length(clust_data.labels)/2)>.5;
    inc_ids = logical([inc_id', ~inc_id']);
    inc_find = find(inc_ids);
    orig_idx = reshape(1:length(clust_data.labels), [2, length(clust_data.labels)/2])';
    final_idx = orig_idx(inc_find);

    ccg_data.pre_cl = ccg_data.pre_cl(final_idx);
    ccg_data.post_cl = ccg_data.post_cl(final_idx);
    clust_data.labels = clust_data.labels(final_idx);
    %% undirected tirin idea
    cat_var = ["cl", "sc", "ct"];
    poses = [   292   198   850   420;     369    72   191   420; 562    71   327   420];
    for k = 1:1
        pre_lab = 'pre_' + cat_var(k);
        post_lab = 'post_' + cat_var(k);
        
        uq_layers = unique(ccg_data.(pre_lab));
        uq_layers = uq_layers(~isnan(uq_layers));
        cur_labels = ccg_data.(cat_var(k) + "_labels");
        pair_labels  = [];
        pair_lay_labels = [];
        cnt = 0;
        for i = 1:length(uq_layers)
            for j = i:length(uq_layers)
                cnt = cnt + 1;
                for kk = 1:clust_data.num_clusters
                    in_layer_pair = (ccg_data.(pre_lab) == uq_layers(i) & ccg_data.(post_lab) == uq_layers(j));% |(ccg_data.(pre_lab) == uq_layers(j) & ccg_data.(post_lab) == uq_layers(i));
                    in_cluster = clust_data.labels == kk;
                    xx = sum(in_layer_pair&in_cluster, 'all');
                    pair_labels(cnt,kk) =  xx;
                end
                pair_lay_labels{cnt} = cur_labels(uq_layers(i)) + "-" + cur_labels(uq_layers(j));
            end
        end
        
        figure('position', poses(k,:));
        x = 1:1:length(pair_lay_labels);
        vq = [];
        xs = [];
        for i = 1:clust_data.num_clusters
            vq(:,i) = interp1(1:length(pair_lay_labels), pair_labels(:,i), x, 'pchip');
            xs(:,i) = x;
        end
        p = plot(xs, vq,'.','linewidth', 2, 'markersize', 20);
        for i = 1:clust_data.num_clusters
            p(i).Color = cmap(i,:);
        end
        set(gca, 'xticklabels', pair_lay_labels, 'xtick', 1:size(pair_labels,1));
        xtickangle(45);
        currx = xlim;
        ylabel("number of pairs");
        xlim([0.5, currx(2)+.5])
        set_axis_defaults();
  end
        save_close_figures(flag.figure_dir + "cluster_layer_counts");

    %% undirected tirin idea prop
    cat_var = ["cl", "sc", "ct"];
    poses = [292   198   850   420;     369    72   191   420; 562    71   327   420];
    for k = 1:1
        pre_lab = 'pre_' + cat_var(k);
        post_lab = 'post_' + cat_var(k);
        
        uq_layers = unique(ccg_data.(pre_lab));
        uq_layers = uq_layers(~isnan(uq_layers));
        cur_labels = ccg_data.(cat_var(k) + "_labels");
        pair_labels  = [];
        pair_lay_labels = [];
        cnt = 0;
        for i = 1:length(uq_layers)
            for j = i:length(uq_layers)
                cnt = cnt + 1;
                for kk = 1:clust_data.num_clusters
                    in_layer_pair = (ccg_data.(pre_lab) == uq_layers(i) & ccg_data.(post_lab) == uq_layers(j));                    in_cluster = clust_data.labels == kk;
                    xx = sum(in_layer_pair&in_cluster, 'all');
                    pair_labels(cnt,kk) =  xx/sum(in_cluster);
                end
                pair_lay_labels{cnt} = cur_labels(uq_layers(i)) + "-" + cur_labels(uq_layers(j));
            end
        end
        
        figure('position', poses(k,:));
        x = 1:1:length(pair_lay_labels);
        vq = [];
        xs = [];
        for i = 1:clust_data.num_clusters
            vq(:,i) = interp1(1:length(pair_lay_labels), pair_labels(:,i), x, 'pchip');
            xs(:,i) = x;
        end
        p = plot(xs, vq,'.','linewidth', 2, 'markersize', 20);
        for i = 1:clust_data.num_clusters
            p(i).Color = cmap(i,:);
        end
        set(gca, 'xticklabels', pair_lay_labels, 'xtick', 1:size(pair_labels,1));
        xtickangle(45);
        currx = xlim;
        ylabel("number of pairs");
        xlim([0.5, currx(2)+.5])
        set_axis_defaults();
    end
    save_close_figures(flag.figure_dir + "cluster_layer_props");
end