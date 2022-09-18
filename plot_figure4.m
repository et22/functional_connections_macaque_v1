%% This script generates plots for figure 4 of the manuscript displaying results from clustering.
%% load cluster output
flag = config();

for group_idx = 1:flag.group_cnt
clearvars -except flag group_idx

load(flag.cluster_output(group_idx), 'clust_data');
load(flag.postproc_output(group_idx), 'ccg_data');

if ~exist(flag.figure_dir, 'dir')
    mkdir(flag.figure_dir)
end

flag.figure_dir = flag.figure_dir + "fig4_";
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


%% new plots
%% silhouette score
figure('position',  [   346   210   586   330]);
plot(clust_data.sil_eva.CriterionValues, 'linewidth', 1.5, 'color', 'k');
hold on;
xline(clust_data.num_clusters, '--', 'linewidth', 1.5, 'color', 'k');
xlim([2,10])
set(gca, 'xtick', [2, 5, 10])
%set(gca, 'ytick', []);
set_axis_defaults()
ylabel("silhouette value")
xlabel("number of clusters")
save_close_figures(flag.figure_dir + 'silscore') 

%% wss score
figure('position',  [   346   210   586   330]);
plot((clust_data.wss(1)-clust_data.wss)/clust_data.wss(1), 'linewidth', 1.5, 'color', 'k');
hold on;
xline(clust_data.num_clusters, '--', 'linewidth', 1.5, 'color', 'k');
xlim([1,10])
set(gca, 'xtick', [1, 5, 10])
%set(gca, 'ytick', []);
set_axis_defaults()
ylabel("silhouette value")
xlabel("number of clusters")
save_close_figures(flag.figure_dir + 'wssscore') 

%% distribution
figure('position',[   97.0000   73.0000  630.6667  554.0000]);
ccg_s = ccg_data.ccgs;
score = clust_data.tsne_mtx;
clust = clust_data.labels(1:length(clust_data.tsne_mtx));

for i = 1:clust_data.num_clusters
    scoren = score(clust == i,:);
    idx = randperm(size(scoren,1));
    
    x = scoren(idx(1:1000),1);
    y = scoren(idx(1:1000),2);
    c = cmap(i,:);
    scatter(x,y,4, 'filled',...
        'MarkerFaceColor', c, 'MarkerEdgeColor', c,...
         'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)
    hold on;
end

xs = get(gca, 'xlim');
xlim([(xs(1) + .1*(xs(2)-xs(1))),xs(2)]);
xs = get(gca, 'xlim');
ys = get(gca, 'ylim');
ylim([ys(1), (ys(2) - .1*(ys(2)-ys(1)))]);
ys = get(gca, 'ylim');

ccg_width = (xs(2)-xs(1))/30;
ccg_height = (ys(2)-ys(1))/20;


for i = 1:clust_data.num_clusters
    num_ccgs = 30;
    ccgn = ccg_s(clust == i, :);
    idxes = ones(num_ccgs,1);
    for j = 1:num_ccgs
        idxes(j) = ceil(rand()*size(ccgn, 1));
    end
    scoren = score(clust == i,:);
    
    ccgs = ccgn(idxes,:);
    ccgs = normalize(ccgs,2); % zscore each row
    rowmin = min(ccgs, [], 2);
    rowmax = max(ccgs, [], 2);
    ccgs = rescale(ccgs, 'inputmin', rowmin, 'inputmax', rowmax);
    posx = scoren(idxes,1)-ccg_width/2;
    posy = scoren(idxes,2)-ccg_height/2;
    
    ccgs = ccgs.*ccg_height + posy;
    ccg_x = [1:1:21]./21.*ccg_width + posx; 
    for ii=1:size(ccgs,1)
        plot(ccg_x(ii,:), ccgs(ii,:), 'color', cmap(i,:));
        hold on;
    end
    hold on;
end

set_axis_defaults
set(gca, 'linewidth', 1)
set(gca, 'Color', 'none') %,'xtick',[], 'ytick',[]
xlabel("t-SNE 1");
ylabel("t-SNE 2");
set(gca, 'Color', 'none', 'XColor', 'none', 'YColor', 'none', 'xtick',[], 'ytick',[])
save_close_figures(flag.figure_dir + 'tsne') 


%% plt ccg templates no axes
for i = 1:clust_data.num_clusters
    figure('position',[   199   430   697/5   170]);
    ccgs_in_cluster = ccg_data.ccgs(clust_data.labels==i,:);
    plot(-flag.max_pt_lag:1:-flag.min_pt_lag,mean(ccgs_in_cluster,1), 'LineWidth', 2, 'Color', cmap(i,:)); % matching convergence
    hold on; 
    %xline(0, 'k--', 'linewidth', 1.5); 
    set(gca, 'linewidth',1.5, 'box', 'off', 'tickdir', 'out', 'xtick', [-flag.max_pt_lag, 0, -flag.min_pt_lag], 'ytick',[]);
    set(gca, 'Color', 'none', 'XColor', 'none', 'YColor', 'none', 'xtick',[], 'ytick',[])
    set_axis_defaults();
    save_close_figures(flag.figure_dir + "ccg_templates_no_axis" + int2str(i)) 
end

%% 
end