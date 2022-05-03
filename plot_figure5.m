%% This script generates plots for figures 5 regarding properties of clusters.
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

%% within vs between proportion comparison
clust_lab = ["sharply synchronous", "broadly synchronous", "forward", "backward"];
for clust_num = 1:4
    cluster = (clust_data.labels == clust_num);
    w_vs_b = ccg_data.pre_cl == ccg_data.post_cl;
    
    n = sum(cluster);
    prop = sum(cluster & w_vs_b)/n;
    exp = .5;
    p = z_prop_test(exp, prop, n);
    disp(clust_lab(clust_num) + " p = " + num2str(p));
end


%% mnr on subset
dists = ccg_data.pair_distance;
within_distance = dists(ccg_data.pre_cl == ccg_data.post_cl);
between_distance = dists(ccg_data.pre_cl ~= ccg_data.post_cl);

subsets{1} = dists > quantile(between_distance, .05) & dists < quantile(within_distance, .95);
%subsets{3} = dists ~= 100009;
%subsets{2} = dists > quantile(between_distance, .1) & dists < quantile(within_distance, .9);

for s = 1:length(subsets)
dists = ccg_data.pair_distance;
within_distance = dists(ccg_data.pre_cl == ccg_data.post_cl);
between_distance = dists(ccg_data.pre_cl ~= ccg_data.post_cl);
subset = subsets{s};
figure; 
hold on;
x_idx = 0;
for clust_num = 1:4
    cluster = categorical(clust_data.labels == clust_num);
    distance = ccg_data.pair_distance;
    w_vs_b = ccg_data.pre_cl == ccg_data.post_cl;

    distance = distance((subset));
    w_vs_b =  w_vs_b(subset);
    cluster = cluster(subset);

    %distance = normalize(distance);
    %w_vs_b = normalize(double(w_vs_b));
    meas = [distance, w_vs_b];

    [B,dev,stats] = mnrfit(meas,cluster);
    if(stats.p(2)<1e-3)
        bar(x_idx+clust_num*2-1, abs(B(2)), 'FaceColor', cmap(clust_num,:), 'EdgeColor', cmap(clust_num,:), 'LineWidth', 1);
    else
        bar(x_idx+clust_num*2-1,abs(B(2)),'FaceColor','w','EdgeColor',cmap(clust_num,:), 'LineWidth', 1);
    end
    errorbar(x_idx+clust_num*2-1,abs(B(2)),stats.se(2), 'Color', 'k', 'linewidth', 1);

    if(stats.p(3)<1e-3)
        bar(x_idx+clust_num*2, abs(B(3)), 'FaceColor', cmap(clust_num,:), 'EdgeColor', cmap(clust_num,:), 'LineWidth', 1);
    else
        bar(x_idx+clust_num*2,abs(B(3)),'FaceColor','w','EdgeColor',cmap(clust_num,:), 'LineWidth', 1);
    end
    errorbar(x_idx+clust_num*2,abs(B(3)),stats.se(3), 'Color', 'k', 'linewidth', 1);
    
    x_idx = x_idx + .5;
    set(gca, 'xtick', []);
    LL = stats.se;% abs(stats.beta - 1.96.*stats.se);
    UL = stats.se;%abs(stats.beta + 1.96.*stats.se);
    disp("subset ana" + clust_lab(clust_num) + " = " + num2str(B(1)) + " + " +...
        num2str(B(2)) + "*distance + " + num2str(B(3)) + "*within_vs_between");
    disp("b conf = " + "[" + num2str(LL(2)) + ", " + num2str(UL(2)) + "]" + "," + ...
        "[" + num2str(LL(3)) + ", " + num2str(UL(3)) + "]");
    disp("p values = " + num2str(stats.p(1)) + ", " + num2str(stats.p(2)) + ", " + num2str(stats.p(3)));
    %disp("p values linhyptest = " + num2str(linhyptest(B, stats.covb, zeros(1,1),[0,1,-1],stats.dfe)));
end
set_axis_defaults();
xl = ylim;
ylim([0, xl(2)]);
end


%% distance matching on subset
dists = ccg_data.pair_distance;
within_distance = dists(ccg_data.pre_cl == ccg_data.post_cl);
between_distance = dists(ccg_data.pre_cl ~= ccg_data.post_cl);
subset = dists > quantile(between_distance, .05) & dists < quantile(within_distance, .95);

cluster = categorical(clust_data.labels == clust_num);
distance = ccg_data.pair_distance;
w_vs_b = ccg_data.pre_cl == ccg_data.post_cl;

distance = distance((subset));
w_vs_b =  w_vs_b(subset);
cluster = cluster(subset);

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
    cluster_prop(1,idx) = sum(simple_to_complex&in_cluster, 'all')/cluster_n;
    cluster_prop(2,idx) = sum(complex_to_simple&in_cluster, 'all')/cluster_n;
    moe(1,idx) = margin_of_error(zvalue, cluster_prop(1,idx), cluster_n);
    moe(2,idx) = margin_of_error(zvalue, cluster_prop(2,idx), cluster_n);
    x = [1,2];
    y = cluster_prop(:,idx);
    err = moe(:,idx);
    errorbar(x+(idx-1-clust_data.num_clusters/2)/20,y,err,'.-', 'color', cmap(idx, :), 'linewidth', 1, 'markersize', 20);
    hold on;
    
    if idx == 3
        X = [sum(simple_to_complex&in_cluster, 'all'), sum(complex_to_simple&in_cluster, 'all')];
        N = [cluster_n, cluster_n];
        [h,p, chi2stat,df] = prop_test(X,N);
        disp("forward s->c vs c ->s prop test: chi = " + num2str(chi2stat)...
            + ", p=" + num2str(p) + ", df=" + num2str(df) + ", n=" + num2str(cluster_n));
    end
end

xlim([0.6,2.4]);
set_axis_defaults();
set(gca, 'xtick', [1, 2], 'xticklabels', ["stc", "cts"])
ys = ylim;
set(gca, 'ytick', [ys(1),ys(1)+(ys(2)-ys(1))/4, ys(1)+(ys(2)-ys(1))/2, ys(1)+3*(ys(2)-ys(1))/4, ys(2)]);
save_close_figures(flag.figure_dir + 'temp51') 

% C-C, S-S vs C-S, S-C test
p = z_prop_test(.5, sum(simple_to_complex + complex_to_simple), length(simple_to_complex));
disp("c-c, s-s vs c-s, s-c: one-proportion z-test: " + num2str(p));

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

disp("within vs between proportions")
disp(cluster_prop);
xlim([0.6,2.4]);
set_axis_defaults();
set(gca, 'xtick', [1, 2], 'xticklabels', ["within", "between"])
ys = ylim;
set(gca, 'ytick', [ys(1),ys(1)+(ys(2)-ys(1))/4, ys(1)+(ys(2)-ys(1))/2, ys(1)+3*(ys(2)-ys(1))/4, ys(2)]);
save_close_figures(flag.figure_dir + 'temp52') 

%% overall chi squared test
layers = unique(ccg_data.(pre_lab));
layers = layers(~isnan(layers));
layer_pairings = zeros(size(ccg_data.(pre_lab)))-1;
for i = 1:length(layers)
    for j = i:length(layers)
        lay_pair = (ccg_data.(pre_lab) == i & ccg_data.(post_lab) == j) |...
        (ccg_data.(pre_lab) == j & ccg_data.(post_lab) == i);
        layer_pairings(lay_pair) = (j) + (i)*length(layers);
    end
end
clusts = clust_data.labels(layer_pairings~=-1);
layer_pairings = layer_pairings(layer_pairings~=-1);
[tbl, chi2, p] = crosstab(clusts, layer_pairings);
dof = (length(unique(clusts))-1)*(length(unique(layer_pairings))-1);
n = length(layer_pairings);
disp("overall layer cross tab: chi2=" + num2str(chi2) + ", p=" + p + ", dof=" + num2str(dof) + ", n=" + int2str(n));

%% pair distance vs cluster
figure('position', [   360   334   320   284]);
groupdata = clust_data.labels;
pds = ccg_data.pair_distance;

b = boxchart(pds, 'GroupByColor', groupdata);
for idx=1:clust_data.num_clusters
    b(idx).BoxFaceColor = cmap(idx,:);
    b(idx).MarkerColor = 'None';
    %b(idx).Orientation = 'Horizontal';
end

ys = [0,1800];
set(gca, 'ylim', ys, 'xtick', [], 'ytick', [ys(1),ys(1)+(ys(2)-ys(1))/4, ys(1)+(ys(2)-ys(1))/2, ys(1)+3*(ys(2)-ys(1))/4, ys(2)]);
set_axis_defaults()
xlabel("pair distance");

save_close_figures(flag.figure_dir + 'temp53') 

clust_labels = ["sharply synchronous", "broadly synchronous", "forward", "backward"];

% clust test
for i = 1:clust_data.num_clusters
    for j = i+1:clust_data.num_clusters
        x1 = pds(clust_data.labels == i);
        x2 = pds(clust_data.labels == j);
        p = ranksum(x1, x2);
        disp("PD Wilcoxon ranksum;" + clust_labels(i) + " vs " + clust_labels(j) + ": p=" + num2str(p));
    end
end

%% r_ori vs cluster
figure('position', [   360   334   320   284]);
groupdata = clust_data.labels;
pds = ccg_data.r_ori;

b = boxchart(pds, 'GroupByColor', groupdata);
for idx=1:clust_data.num_clusters
    b(idx).BoxFaceColor = cmap(idx,:);
    b(idx).MarkerColor = 'None';
    %b(idx).Orientation = 'Horizontal';
end
%ys = [0,1800];
%set(gca, 'xlim', ys, 'ytick', [], 'xtick', [ys(1),ys(1)+(ys(2)-ys(1))/4, ys(1)+(ys(2)-ys(1))/2, ys(1)+3*(ys(2)-ys(1))/4, ys(2)]);
set_axis_defaults()
xlabel("r_ori");

save_close_figures(flag.figure_dir + 'temp54') 


% clust test
for i = 1:clust_data.num_clusters
    for j = i+1:clust_data.num_clusters
        x1 = pds(clust_data.labels == i);
        x2 = pds(clust_data.labels == j);
        p = ranksum(x1, x2);
        disp("Rori Wilcoxon ranksum;" + clust_labels(i) + " vs " + clust_labels(j) + ": p=" + num2str(p));
    end
end
end

function y = margin_of_error(z_value, prop, n)
    y =  z_value*sqrt((prop*(1-prop))/n);
end