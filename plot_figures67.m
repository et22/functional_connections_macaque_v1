%% This script generates plots for figures 6-7 regarding the distribution of clusters across layers.
%% load cluster output
flag = config();

for group_idx = 1:flag.group_cnt
clearvars -except flag group_idx

load(flag.cluster_output(group_idx), 'clust_data');
load(flag.postproc_output(group_idx), 'ccg_data');

if ~exist(flag.figure_dir, 'dir')
    mkdir(flag.figure_dir)
end

flag.figure_dir = flag.figure_dir + "fig67_";

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

% generate random half
rng(2); % seed rng for replicability
inc_id = rand(1,length(clust_data.labels)/2)>.5;
inc_ids = logical([inc_id', ~inc_id']);
inc_find = find(inc_ids);
orig_idx = reshape(1:length(clust_data.labels), [2, length(clust_data.labels)/2])';
final_idx = orig_idx(inc_find);

%% new code 
%overall params
cl_labels = ["2/3", "4a/b", "4c\alpha", "4c\beta", "5", "6", "WM"];
sc_labels = ccg_data.sc_labels;
moe_crit = .95;
chi_crit = .05/24;
lab = 'cl';
% figure 6
%spec params
pre_layers = ["4c\alpha"];
post_layers = ["4c\alpha"];
figpath="4ca + 4ca";
half_sub = true;
plot_subplot(pre_layers, post_layers, figpath, moe_crit, chi_crit, cl_labels, ccg_data, cmap, clust_data, half_sub, final_idx, lab)

%spec params
pre_layers = ["4c\beta"];
post_layers = ["4c\beta"];
figpath="4cb+4cb";
half_sub = true;
plot_subplot(pre_layers, post_layers, figpath, moe_crit, chi_crit, cl_labels, ccg_data, cmap, clust_data, half_sub, final_idx, lab)

%spec params
pre_layers = ["4c\beta", "4c\alpha"];
post_layers = ["4c\alpha", "4c\beta",];
figpath="4ca + 4cb";
half_sub = true;
plot_subplot(pre_layers, post_layers, figpath, moe_crit, chi_crit, cl_labels, ccg_data, cmap, clust_data, half_sub, final_idx, lab)

%spec params
pre_layers = ["6"];
post_layers = ["4c\alpha"];
figpath="6 to 4ca";
half_sub = false;
plot_subplot(pre_layers, post_layers, figpath, moe_crit, chi_crit, cl_labels, ccg_data, cmap, clust_data, half_sub, final_idx, lab)

%spec params
pre_layers = ["6"];
post_layers = ["4c\beta"];
figpath="6 to 4cb";
half_sub = false;
plot_subplot(pre_layers, post_layers, figpath, moe_crit, chi_crit, cl_labels, ccg_data, cmap, clust_data, half_sub, final_idx, lab)

%spec params
pre_layers = ["6"];
post_layers = ["6"];
figpath="6 + 6";
half_sub = true;
plot_subplot(pre_layers, post_layers, figpath, moe_crit, chi_crit, cl_labels, ccg_data, cmap, clust_data, half_sub, final_idx, lab)


% figure 7
%spec params
pre_layers = ["6"];
post_layers = ["2/3"];
figpath="6 to 23";
half_sub = false;
plot_subplot(pre_layers, post_layers, figpath, moe_crit, chi_crit, cl_labels, ccg_data, cmap, clust_data, half_sub, final_idx, lab)

%spec params
pre_layers = ["5"];
post_layers = ["2/3"];
figpath="5 to 23";
half_sub = false;
plot_subplot(pre_layers, post_layers, figpath, moe_crit, chi_crit, cl_labels, ccg_data, cmap, clust_data, half_sub, final_idx, lab)

%spec params
pre_layers = ["4c\beta"];
post_layers = ["2/3"];
figpath="4cb to 23";
half_sub = false;
plot_subplot(pre_layers, post_layers, figpath, moe_crit, chi_crit, cl_labels, ccg_data, cmap, clust_data, half_sub, final_idx, lab)

%spec params
pre_layers = ["4c\alpha"];
post_layers = ["2/3"];
figpath="4ca to 23";
half_sub = false;
plot_subplot(pre_layers, post_layers, figpath, moe_crit, chi_crit, cl_labels, ccg_data, cmap, clust_data, half_sub, final_idx, lab)

%spec params
pre_layers = ["4a/b"];
post_layers = ["2/3"];
figpath="4ab to 23";
half_sub = false;
plot_subplot(pre_layers, post_layers, figpath, moe_crit, chi_crit, cl_labels, ccg_data, cmap, clust_data, half_sub, final_idx, lab)

%spec params
pre_layers = ["2/3"];
post_layers = ["2/3"];
figpath="23 + 23";
half_sub = true;
plot_subplot(pre_layers, post_layers, figpath, moe_crit, chi_crit, cl_labels, ccg_data, cmap, clust_data, half_sub, final_idx, lab)

end


function plot_subplot(pre_layers, post_layers, figpath, moe_crit, chi_crit, cl_labels, ccg_data, cmap, clust_data, half_sub, final_idx, lab)
clust_labs = ["ssync", "bsync", "fasync", "rasync"];
disp("=================================================================");
disp("overall cross tab")
pre_lab = append('pre_', lab);
post_lab = append('post_', lab);
conf_threshold = moe_crit;

if conf_threshold == .95
    zvalue = 1.96;
else
    error("undefined confidence value")
end

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
figure('position', [  360.0000  445.0000  240.3333  173.0000]);
for idx = 1:clust_data.num_clusters
    in_cluster = clust_data.labels == idx;
    
    in_cluster = in_cluster(final_idx);
    cluster_n = sum(in_cluster, 'all');

    cluster_prob(idx) = sum(inlayerpair&in_cluster, 'all')/cluster_n;
    moe(idx) = margin_of_error(zvalue, cluster_prob(idx), cluster_n);
    x = idx;
    y = cluster_prob(idx);
    err = moe(idx);
    errorbar(x,y*100,err*100,'.-', 'color', cmap(idx, :), 'linewidth', 1, 'markersize', 20);
    hold on;
end

xlim([0.6,clust_data.num_clusters + .4]);
set_axis_defaults();
set(gca, 'xtick', [])
ys = ylim;
ylim([0,ys(2)]);

ys = ylim;
ticklabs =  string(round([ys(1), ys(1)+(ys(2)-ys(1))/2, ys(2)],2));
ticklabs = equalize_label_lengths(ticklabs);

set(gca, 'ytick', round([ys(1), ys(1)+(ys(2)-ys(1))/2, ys(2)],2), 'yticklabels', ticklabs);


% sig testing
% add overall n to plots for that subdivision
text(.7,.7,"n = " + int2str(sum(inlayerpair)), 'units', 'normalized', 'fontsize', 11);

clust_data.labels = clust_data.labels(final_idx);

%overall chi2 test
in_pair = inlayerpair;
[tbl, chi2, p] = crosstab(clust_data.labels, in_pair);
dof = (length(unique(clust_data.labels))-1)*1;
n = sum(inlayerpair);
disp(figpath + " cross tab: chi2=" + num2str(chi2) + ", p=" + p + ", dof=" + num2str(dof) + ", n=" + int2str(n));

if p<chi_crit
    text(1,1,'*', 'units', 'normalized', 'fontsize', 20);
end

%pairwise chi2 test
y = cluster_prob;
[m, i] = max(y);
ys = ylim;
yrange = ys(2)-ys(1);

disp("pairwise comparisons with " + clust_labs(i));
for idx = 1:clust_data.num_clusters
    sub_idx = clust_data.labels == i | clust_data.labels == idx;
    [tbl, chi2, p] = crosstab(clust_data.labels(sub_idx), in_pair(sub_idx));
    if p<chi_crit
        text(idx+.2,cluster_prob(idx)*100 + .02*yrange,'*', 'units', 'data','horizontalalignment', 'center','fontsize', 14);
    end
    if idx ~= i
        disp(figpath + " cluster " + clust_labs(idx) + " cross tab: chi2=" + num2str(chi2) + ", p=" + p);
    end
end

ylabel(figpath + " (%)");


save_close_figures("figures/fig67_" + figpath) 

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