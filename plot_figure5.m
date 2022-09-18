%% This script generates plots for figures 5 regarding properties of clusters.
%% load cluster output
addpath(genpath(pwd));

flag = config();

for group_idx = 1:flag.group_cnt
clearvars -except flag group_idx

load(flag.cluster_output(group_idx), 'clust_data');
load(flag.postproc_output(group_idx), 'ccg_data');

if ~exist(flag.figure_dir, 'dir')
    mkdir(flag.figure_dir)
end

flag.figure_dir = flag.figure_dir + "fig5_";

%% subsetting ccg_data based on significance
ccg_data_all = ccg_data.ccg;
[ccg_data, sig_idx] = get_significant_ccgs(ccg_data_all, flag);

ccg_data.cl_labels = ["2/3", "4a/b", "4c\alpha", "4c\beta", "5", "6", "WM"];
ccg_data.sc_labels = ["comp.", "simp."];
ccg_data.ct_labels = ["AS", "FS", "RM", "RL"];

%% set up plotting
cmap = [14,121,178;243,146,55;191,19,99;244,159,188;6,214,160]./255;%[linspace(247, 198, clust_data.num_clusters)', linspace(219,121, clust_data.num_clusters)', linspace(187,88,clust_data.num_clusters)']./255;
cmap = [cmap; 1-cmap];
cmap=cmap(1:clust_data.num_clusters,:);
[empty_label{1:clust_data.num_clusters}] = deal('');
clust_data.labels = sort_clusters_by_lag(clust_data.num_clusters, ccg_data.ccgs, clust_data.labels);



%% generate final_idx, random direction for CCGs to use for figures
% generate random half
rng(2); % seed rng for replicability
inc_id = rand(1,length(clust_data.labels)/2)>.5;
inc_ids = logical([inc_id', ~inc_id']);
inc_find = find(inc_ids);
orig_idx = reshape(1:length(clust_data.labels), [2, length(clust_data.labels)/2])';
final_idx = orig_idx(inc_find);


%% pair distance vs cluster (random half)
disp("===========pair distance vs cluster==================")
figure('position', [   360   334   320   284]);

groupdata = clust_data.labels(final_idx);
pds = ccg_data.pair_distance(final_idx);

b = boxchart(pds, 'GroupByColor', groupdata);
for idx=1:clust_data.num_clusters
    b(idx).BoxFaceColor = cmap(idx,:);
    b(idx).MarkerColor = 'None';
end

ys = [0,1800];
set(gca, 'ylim', ys, 'xtick', [], 'ytick', [ys(1),ys(1)+(ys(2)-ys(1))/4, ys(1)+(ys(2)-ys(1))/2, ys(1)+3*(ys(2)-ys(1))/4, ys(2)]);
set_axis_defaults()
xlabel("pair distance");

clust_labels = ["sharply synchronous", "broadly synchronous", "forward", "backward"];

% clust test
for i = 1:clust_data.num_clusters
    for j = i+1:clust_data.num_clusters
        x1 = pds(groupdata == i);
        x2 = pds(groupdata == j);
        p = ranksum(x1, x2);
        disp("PD Wilcoxon ranksum;" + clust_labels(i) + " vs " + clust_labels(j) + ": p=" + num2str(p));
    end
end

save_close_figures(flag.figure_dir + 'pd_v_clust') 

%% r_ori vs cluster half
disp("===========r ori vs cluster==================")
figure('position', [   360   334   320   284]);

groupdata = clust_data.labels(final_idx);
pds = ccg_data.r_ori(final_idx);

b = boxchart(pds, 'GroupByColor', groupdata);
for idx=1:clust_data.num_clusters
    b(idx).BoxFaceColor = cmap(idx,:);
    b(idx).MarkerColor = 'None';
end
set_axis_defaults()
xlabel("r_ori");

% clust test
for i = 1:clust_data.num_clusters
    for j = i+1:clust_data.num_clusters
        x1 = pds(groupdata == i);
        x2 = pds(groupdata == j);
        p = ranksum(x1, x2);
        disp("Rori Wilcoxon ranksum;" + clust_labels(i) + " vs " + clust_labels(j) + ": p=" + num2str(p));
    end
end

save_close_figures(flag.figure_dir + 'rori_v_clust') 

%% layer pair vs cluster cross tab
disp("===========layer pair vs cluster==================")
pre_cl = ccg_data.pre_cl(final_idx);
post_cl = ccg_data.post_cl(final_idx);
clusts = clust_data.labels(final_idx);
layers = unique(pre_cl);
cnt = 0;
for k = 1:length(pre_cl)
    if ~isnan(pre_cl(k)) && ~isnan(post_cl(k))
        cnt = cnt + 1;
        i1 = find(layers == pre_cl(k));
        i2 = find(layers == post_cl(k));
        lay_pair(cnt) = max(i1,i2)+min(i1,i2)*(length(layers)+1);
        clus(cnt) = clusts(k);
    end
end
[tbl,chi2,p] = crosstab(lay_pair, clus);
disp("Cross-tab layer-pair vs cluster overall: p=" + num2str(p));

%% logistic regression on subset of data
disp("===========logistic regression analysis==================")
clust_lab = ["sharply synchronous", "broadly synchronous", "asynchronous"];
rng(1);
dists = ccg_data.pair_distance(final_idx);
within_distance = dists(ccg_data.pre_cl(final_idx) == ccg_data.post_cl(final_idx));
between_distance = dists(ccg_data.pre_cl(final_idx) ~= ccg_data.post_cl(final_idx));

subset = dists > quantile(between_distance, .05) & dists < quantile(within_distance, .95);

x_idx = 0;
for clust_num = 1:3
    if clust_num < 3
        cluster = categorical(clust_data.labels(final_idx) == clust_num);
    else
        cluster = categorical(clust_data.labels(final_idx) == 3 | clust_data.labels(final_idx) == 4);
    end
    distance = ccg_data.pair_distance(final_idx);
    w_vs_b = ccg_data.pre_cl(final_idx) == ccg_data.post_cl(final_idx);

    distance = distance((subset));
    w_vs_b =  w_vs_b(subset);
    cluster = cluster(subset);

    meas = [distance, w_vs_b];

    [B,dev,stats] = mnrfit(meas,cluster);
    
    LL = round(stats.se, 4);
    B = round(B,4);
    
    disp("" + clust_lab(clust_num) + " = " + num2str(B(1)) + " + " +...
        num2str(B(2)) + "*distance + " + num2str(B(3)) + "*within_vs_between");
    disp("b conf = " + "[" + num2str(LL(2)) + "]" + "," + ...
        "[" + num2str(LL(3)) + "]");
    disp("p values = " + num2str(stats.p(1)) + ", " + num2str(stats.p(2)) + ", " + num2str(stats.p(3)));
end


%% within vs between proportion comparison
disp("===========within vs between==================")
clust_lab = ["sharply synchronous", "broadly synchronous", "forward", "backward"];
for clust_num = 1:4
    cluster = (clust_data.labels(final_idx) == clust_num);
    w_vs_b = ccg_data.pre_cl(final_idx) == ccg_data.post_cl(final_idx);
    
    n = sum(cluster);
    prop = sum(cluster & w_vs_b)/n;
    exp = .5;
    p = z_prop_test(exp, prop, n);
    disp(clust_lab(clust_num) + "prop. = " + num2str(prop) + " p = " + num2str(p));
end


%% plot within vs between (half)
pre_lab = 'pre_cl';
post_lab = 'post_cl';
conf_threshold = .95;

if conf_threshold == .95
    zvalue = 1.96;
else
    error("undefined confidence value")
end

for idx = 1:clust_data.num_clusters
    within = ccg_data.(pre_lab) == ccg_data.(post_lab);
    between = ccg_data.(pre_lab) ~= ccg_data.(post_lab);
    in_cluster = clust_data.labels == idx;

    % half index
    within = within(final_idx);
    between = between(final_idx);
    in_cluster = in_cluster(final_idx);

    cluster_n = sum(in_cluster, 'all');
    cluster_prop(1,idx) = sum(within&in_cluster, 'all')/cluster_n;
    cluster_prop(2,idx) = sum(between&in_cluster, 'all')/cluster_n;
    moe(1,idx) = margin_of_error(zvalue, cluster_prop(1,idx), cluster_n);
    moe(2,idx) = margin_of_error(zvalue, cluster_prop(2,idx), cluster_n);
    x = [1,2];
    y = cluster_prop(:,idx);
    err = moe(:,idx);
end


disp("within vs between proportions")
disp(cluster_prop);

%save_close_figures(flag.figure_dir + 'withinbetweenhalf') 
figure('position', [   360   334   289   284]);
hold on;
for j = 1:size(cluster_prop, 2)
    b1 = bar(j-.2, cluster_prop(1,j),.4, 'FaceColor',cmap(j,:), 'FaceAlpha', .2, 'EdgeColor', cmap(j,:), 'linewidth', 1);
    b2 = bar(j+.2, cluster_prop(2,j),.4, 'FaceColor',cmap(j,:), 'FaceAlpha', .2, 'EdgeColor', cmap(j,:), 'linewidth', 1);
    hatchfill2(b2,'single','HatchAngle',45,'hatchcolor',cmap(j,:));

    %for k = 1:5
    %    scatter(j, cprop_ses(1,j,k)-cprop_ses(2,j,k), 30, 'k', 'filled');
    %end
end
for i = 1:length(p)
    if p(i)<1e-5
        text(i, cluster_prop(1,i)-cluster_prop(2,i), '*', 'horizontalalignment', 'center');
    end
end
set_axis_defaults();
set(gca, 'xtick', [1, 2, 3, 4], 'xticklabels', ["1", "2", "3", "4"])
xlim([0.4, 4.6]);

save_close_figures(flag.figure_dir + 'within_vs_between') 

%% simple-simple vs simple-complex vs complex-complex z-test
disp("===========simples vs complex z-test==================")
simp_simp = sum(ccg_data.pre_sc(final_idx) == 2 & ccg_data.post_sc(final_idx) == 2);
comp_comp = sum(ccg_data.pre_sc(final_idx) == 1 & ccg_data.post_sc(final_idx) == 1);
simp_comp = sum((ccg_data.pre_sc(final_idx) == 1 & ccg_data.post_sc(final_idx) == 2) | ...
              (ccg_data.pre_sc(final_idx) == 2 & ccg_data.post_sc(final_idx) == 1));
all_n = simp_simp + comp_comp + simp_comp;
disp("proportions s-s: " + simp_simp/all_n + ", c-c: " + comp_comp/all_n + ", s-c: " + simp_comp/all_n);

p = z_prop_test(.5, comp_comp/all_n, all_n);
disp("one prop. z test: p=" + p);

%% simple vs complex (half)
disp("===========simple-complex paired tests==================")
%overall params
cl_labels = ["2/3", "4a/b", "4c\alpha", "4c\beta", "5", "6", "WM"];
sc_labels = ["comp.", "simp."];

%spec params
pre_layers = ["simp."];
post_layers = ["comp."];
figpath="S to C";
half_sub = false;
lab = "sc";
plot_subplot(pre_layers, post_layers, figpath, sc_labels, ccg_data, cmap, clust_data, half_sub, final_idx, lab);
save_close_figures(flag.figure_dir + 'simple_vs_complex') 
end

function plot_subplot(pre_layers, post_layers, figpath, cl_labels, ccg_data, cmap, clust_data, half_sub, final_idx, lab)
pre_lab = append('pre_', lab);
post_lab = append('post_', lab);

inlayerpair = zeros(size(ccg_data.(pre_lab)));
for i = 1:length(pre_layers) 
    pre_l = find(cl_labels == pre_layers(i));
    post_l = find(cl_labels == post_layers(i));
    inlayerpair = inlayerpair | (ccg_data.(pre_lab) == pre_l & ccg_data.(post_lab) == post_l);
end

%% 
if half_sub
    inlayerpair = inlayerpair(final_idx);
else
    inlayerf = find(inlayerpair);
    notin = setdiff(inlayerf, final_idx);
    for i = 1:length(notin)
        if (mod(notin(i),2) == 0)
            final_idx(final_idx==(notin(i)-1)) = notin(i);
        else
            final_idx(final_idx==(notin(i)+1)) = notin(i);
        end
    end
    inlayerpair = inlayerpair(final_idx);
end

%% plotting code
for idx = 1:clust_data.num_clusters
    in_cluster = (clust_data.labels(final_idx) == idx) & inlayerpair;

    all_n = sum(inlayerpair);

    cluster_prob(idx) = sum(in_cluster, 'all')/all_n;
    cluster_mem(idx) = sum(in_cluster, 'all');
    x = idx;
    y = cluster_prob(idx);
end

figure('position', [   360   334   289   284]);
hold on;
for j = 1:length(cluster_prob)
    b1 = bar(j, cluster_prob(j),.75, 'FaceColor',cmap(j,:), 'FaceAlpha', .2, 'EdgeColor', cmap(j,:), 'linewidth', 1);
end
set_axis_defaults();
set(gca, 'xtick', [1, 2, 3, 4], 'xticklabels', ["1", "2", "3", "4"])

% sig testing
% add overall n to plots for that subdivision
text(.7,.7,"n = " + int2str(sum(inlayerpair)), 'units', 'normalized', 'fontsize', 11);

clust_data.labels = clust_data.labels(final_idx);
clust_labs = ["ssync", "bsync", "fasync", "rasync"];
disp("expectation is uniform:")
for i = 1:clust_data.num_clusters
    for j = i+1:clust_data.num_clusters
        bins = (0:1)';
        obsCounts = [cluster_mem(i); cluster_mem(j)];
        numi = sum(clust_data.labels == i);
        numj = sum(clust_data.labels == j);
        expCounts = [(sum(obsCounts)); (sum(obsCounts))];
        [h, p, stats] = chi2gof(bins, 'Ctrs', bins, 'Frequency', obsCounts', 'Expected', expCounts);
        disp(clust_labs(i) + " vs " + clust_labs(j) + " p = " + num2str(p));
    end
end

disp("expectation is proportional to size of class:")
for i = 1:clust_data.num_clusters
    for j = i+1:clust_data.num_clusters
        bins = (0:1)';
        obsCounts = [cluster_mem(i); cluster_mem(j)];
        numi = sum(clust_data.labels == i);
        numj = sum(clust_data.labels == j);
        itoj = numi/(numi+numj);
        jtoi = numj/(numi+numj);
        expCounts = [itoj*(sum(obsCounts)); jtoi*(sum(obsCounts))];
        [h, p, stats] = chi2gof(bins, 'Ctrs', bins, 'Frequency', obsCounts', 'Expected', expCounts);
        disp(clust_labs(i) + " vs " + clust_labs(j) + " p = " + num2str(p));
    end
end

end

function y = margin_of_error(z_value, prop, n)
    y =  z_value*sqrt((prop*(1-prop))/n);
end

function ticklabs = equalize_label_lengths(ticklabs)
for lab = 1:length(ticklabs)
    while length(char(ticklabs(lab)))<4
        if contains(ticklabs(lab),".")
            ticklabs(lab) = ticklabs(lab) + "0";
        else
            ticklabs(lab) = ticklabs(lab) + ".";
        end
    end
    while length(char(ticklabs(lab)))>4
        charvers = char(ticklabs(lab));
        charvers = charvers(1:end-1);
        ticklabs(lab) = string(charvers);
    end
end
end