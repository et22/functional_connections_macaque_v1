%% load cluster output
flag = config();

for group_idx = 1:flag.group_cnt
clearvars -except flag group_idx
if ~flag.augment
    load(flag.cluster_output(group_idx), 'clust_data');
else
    load(flag.augment_output, 'clust_data');
end
load(flag.postproc_output(group_idx), 'ccg_data');

if ~exist(flag.figure_dir, 'dir')
    mkdir(flag.figure_dir)
end

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

%% plot simple vs complex 
pre_lab = 'pre_sc';
post_lab = 'post_sc';
conf_threshold = .95;

if conf_threshold == .95
    zvalue = 1.96;
else
    error("undefined confidence value")
end

figure('position', [   360   334   289   284]);
for idx = 1:clust_data.num_clusters
    complex_to_simple = ccg_data.(pre_lab) == 1 & ccg_data.(post_lab) == 2;
    simple_to_complex = ccg_data.(pre_lab) == 2 & ccg_data.(post_lab) == 1;
    in_cluster = clust_data.labels == idx;
    cluster_n = sum(in_cluster, 'all');
    cluster_prop(1,idx) = sum(complex_to_simple&in_cluster, 'all')/cluster_n;
    cluster_prop(2,idx) = sum(simple_to_complex&in_cluster, 'all')/cluster_n;
    moe(1,idx) = margin_of_error(zvalue, cluster_prop(1,idx), cluster_n);
    moe(2,idx) = margin_of_error(zvalue, cluster_prop(2,idx), cluster_n);
    x = [1,2];
    y = cluster_prop(:,idx);
    err = moe(:,idx);
    errorbar(x+(idx-1-clust_data.num_clusters/2)/20,y,err,'.-', 'color', cmap(idx, :), 'linewidth', 1, 'markersize', 20);
    hold on;
end

xlim([0.6,2.4]);
set_axis_defaults();
set(gca, 'xtick', [1, 2], 'xticklabels', ["cts", "stc"])
ys = ylim;
set(gca, 'ytick', [ys(1),ys(1)+(ys(2)-ys(1))/4, ys(1)+(ys(2)-ys(1))/2, ys(1)+3*(ys(2)-ys(1))/4, ys(2)]);
save_close_figures('temp51') 



%% plot within vs between 
pre_lab = 'pre_cl';
post_lab = 'post_cl';
conf_threshold = .95;

if conf_threshold == .95
    zvalue = 1.96;
else
    error("undefined confidence value")
end

figure('position', [   360   334   289   284]);
for idx = 1:clust_data.num_clusters
    within = ccg_data.(pre_lab) == ccg_data.(post_lab);
    between = ccg_data.(pre_lab) ~= ccg_data.(post_lab);
    in_cluster = clust_data.labels == idx;
    cluster_n = sum(in_cluster, 'all');
    cluster_prop(1,idx) = sum(within&in_cluster, 'all')/cluster_n;
    cluster_prop(2,idx) = sum(between&in_cluster, 'all')/cluster_n;
    moe(1,idx) = margin_of_error(zvalue, cluster_prop(1,idx), cluster_n);
    moe(2,idx) = margin_of_error(zvalue, cluster_prop(2,idx), cluster_n);
    x = [1,2];
    y = cluster_prop(:,idx);
    err = moe(:,idx);
    errorbar(x+(idx-1-clust_data.num_clusters/2)/20,y,err,'.-', 'color', cmap(idx, :), 'linewidth', 1, 'markersize', 20);
    hold on;
end

xlim([0.6,2.4]);
set_axis_defaults();
set(gca, 'xtick', [1, 2], 'xticklabels', ["within", "between"])
ys = ylim;
set(gca, 'ytick', [ys(1),ys(1)+(ys(2)-ys(1))/4, ys(1)+(ys(2)-ys(1))/2, ys(1)+3*(ys(2)-ys(1))/4, ys(2)]);
save_close_figures('temp52') 

%% pair distance vs cluster
figure('position', [   360   334   289   284]);
groupdata = clust_data.labels;
pds = ccg_data.pair_distance;

b = boxchart(pds, 'GroupByColor', groupdata);
for idx=1:clust_data.num_clusters
    b(idx).BoxFaceColor = cmap(idx,:);
    b(idx).MarkerColor = 'None';
    b(idx).Orientation = 'Horizontal';
end
ys = [0,1800];
set(gca, 'xlim', ys, 'ytick', [], 'xtick', [ys(1),ys(1)+(ys(2)-ys(1))/4, ys(1)+(ys(2)-ys(1))/2, ys(1)+3*(ys(2)-ys(1))/4, ys(2)]);
set_axis_defaults()
xlabel("pair distance");

save_close_figures('temp53') 

%% r_ori vs cluster
figure('position', [   360   334   289   284]);
groupdata = clust_data.labels;
pds = ccg_data.pair_distance;

b = boxchart(pds, 'GroupByColor', groupdata);
for idx=1:clust_data.num_clusters
    b(idx).BoxFaceColor = cmap(idx,:);
    b(idx).MarkerColor = 'None';
    b(idx).Orientation = 'Horizontal';
end
ys = [0,1800];
set(gca, 'xlim', ys, 'ytick', [], 'xtick', [ys(1),ys(1)+(ys(2)-ys(1))/4, ys(1)+(ys(2)-ys(1))/2, ys(1)+3*(ys(2)-ys(1))/4, ys(2)]);
set_axis_defaults()
xlabel("r_ori");

save_close_figures('temp54') 

end

function y = margin_of_error(z_value, prop, n)
    y =  z_value*sqrt((prop*(1-prop))/n);
end