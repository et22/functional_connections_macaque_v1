%% load cluster output
flag = config();
load(flag.cluster_output, 'clust_data');
load(flag.postproc_output, 'ccg_data');


%% subsetting ccg_data based on significance
ccg_data = ccg_data.ccg;
sig_idx = (ccg_data.noise_std2>flag.sig_min_std ) & ...
        (ccg_data.peaks>(flag.sig_num_stds*ccg_data.noise_std2 + ccg_data.noise_mean2)) & ...
        (abs(ccg_data.peak_lag) <= flag.sig_max_lag);

fields = fieldnames(ccg_data);
for i = 1:length(fields)
    if ~strcmp(fields{i},'cluster') && ~strcmp(fields{i},'config') && ~strcmp(fields{i},'ccg_control')
        ccg_data.(fields{i}) = ccg_data.(fields{i})(sig_idx,:);
    end
end

ccg_data.cl_labels = ["2/3", "4a/b", "4c\alpha", "4c\beta", "5", "6", "WM"];
ccg_data.sc_labels = ["Complex", "Simple"];
ccg_data.ct_labels = ["AS", "FS", "RM", "RL"];
%% set up plotting
cmap = [linspace(247, 198, clust_data.num_clusters)', linspace(219,121, clust_data.num_clusters)', linspace(187,88,clust_data.num_clusters)']./255;
clust_data.labels = sort_clusters_by_lag(clust_data.num_clusters, ccg_data.ccgs, clust_data.labels);

%% plot cluster validation metrics
figure('position',  [   346   210   586   330]);
subplot(1,2,1); % sil method
plot(clust_data.sil_eva.CriterionValues, 'linewidth', 1.5, 'color', 'k');
hold on;
xline(clust_data.num_clusters, '--', 'linewidth', 1.5, 'color', 'k');
xlim([2,10])
set(gca, 'xtick', [2, 5, 10])
%set(gca, 'ytick', []);
set_axis_defaults()
ylabel("silhouette value")
xlabel("number of clusters")

% calculating w_k
for i = 1:10
    wss(i) = w_k(clust_data.tsne_mtx, clust_data.clusters(:,i));
end
expl_var = (wss(1)-wss)/wss(1); 
subplot(1,2,2); % elbow

plot(expl_var*100, 'linewidth', 1.5, 'color', 'k');
hold on;
xline(clust_data.num_clusters, '--', 'linewidth', 1.5, 'color', 'k');
xlim([1,10])
set(gca, 'xtick', [1,4,10])
%set(gca, 'ytick', []);
set_axis_defaults()
ylabel("explained variance (%)")
xlabel("number of clusters")
save_close_figures(flag.figure_dir + "silhouette") 

%% plot cluster in 3D tSNE space
figure('position', [445   113   457   387]);
for i = 1:clust_data.num_clusters
    curr_tsne_score = clust_data.tsne_mtx(clust_data.labels==i,:);
    scatter3(curr_tsne_score(:,3), curr_tsne_score(:,1),curr_tsne_score(:,2),4, 'filled', 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor', cmap(i,:))
    hold on;
end
set_axis_defaults
xlabel("t-SNE 3");
ylabel("t-SNE 1");
zlabel("t-SNE 2");
set(gca, 'view', [-70.3000    9.9805]);
save_close_figures(flag.figure_dir + "3dscatter_tsne") 

%% plot cluster in 2D tSNE space w/out axes
%figure('position', [445   113   457   387]);
figure('position',[97    32   645   595]);
for i = 1:clust_data.num_clusters
    curr_tsne_score = clust_data.tsne_mtx(clust_data.labels==i,:);
    scatter(curr_tsne_score(:,1),curr_tsne_score(:,2),4, 'filled', 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor', cmap(i,:))
    hold on;
end
set_axis_defaults
set(gca, 'Color', 'none', 'XColor', 'none', 'YColor', 'none', 'xtick',[], 'ytick',[])
save_close_figures(flag.figure_dir + "2dscatter_tsne") 

%% plot attributes x clusters
plot_attributes_x_cluster(ccg_data,clust_data, cmap);
save_close_figures(flag.figure_dir + "ccg_attributes") 

%% plt ccg templates with axes
figure('position',[   199   430   697   170]);
for i = 1:clust_data.num_clusters
    ccgs_in_cluster = ccg_data.ccgs(clust_data.labels==i,:);
    subplot(1,clust_data.num_clusters,i);
    plot(-flag.max_pt_lag:1:-flag.min_pt_lag,mean(ccgs_in_cluster,1), 'LineWidth', 2, 'Color', cmap(i,:)); % matching convergence
    hold on; xline(0, 'k--', 'linewidth', 1.5); 
    set(gca, 'linewidth',1.5, 'box', 'off', 'tickdir', 'out', 'xtick', [-flag.max_pt_lag, 0, -flag.min_pt_lag], 'ytick',[]);
    xlabel('\tau');
    ylabel('');
    set(gca, 'Color', 'none', 'YColor', 'none', 'ytick',[])
    text(.9,1.05,"n = "+int2str(sum(clust_data.labels==i)), 'horizontalalignment', 'center','units', 'normalized', 'fontsize', 10);

    set_axis_defaults();
end
save_close_figures(flag.figure_dir + "ccg_templates") 

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

%% plot example ccgs
figure('position',[187    30   699   609]);
rng(flag.rng_seed+2) 
for i = 1:clust_data.num_clusters
    num_cluster_mems = sum(clust_data.labels==i);
    ccgs_in_cluster = ccg_data.ccgs(clust_data.labels==i,:);
    for j = 1:3
        labels1 = ceil(rand()*num_cluster_mems); 
        subplot(3,clust_data.num_clusters,i + (j-1)*clust_data.num_clusters);
        hold on; xline(0, 'k--', 'linewidth', 1.5); 
        plot(-flag.max_pt_lag:1:-flag.min_pt_lag,mean(ccgs_in_cluster(labels1,:),1), 'LineWidth', 2, 'Color', cmap(i,:)); % matching convergence
        set(gca, 'linewidth',1.5, 'box', 'off', 'tickdir', 'out', 'xtick', [-flag.max_pt_lag, 0, -flag.min_pt_lag], 'ytick',[]);
        xlabel('\tau');
        ylabel('');
        set(gca, 'Color', 'none', 'YColor', 'none', 'ytick',[])
        set_axis_defaults()
    end
end
save_close_figures(flag.figure_dir + "ccg_examples") 

%% crosstabs within/between
wb_layers = ccg_data.pre_cl == ccg_data.post_cl;
[conttbl,chi2,p,labels] = crosstab(wb_layers, clust_data.labels);
cluster_cnts = sum(conttbl,1);
wb_cnts = sum(conttbl,2);
conttbl_prop = conttbl./cluster_cnts;
disp("x-tab within vs between layers: " + "chi^2 = " + chi2 + ", p=" + p);

figure('position', [360   198   560   300]);
h = heatmap(round(conttbl,2), 'Colormap', bone);
ax = gca();
h.XDisplayLabels = nan(size(ax.XDisplayData));
h.YDisplayLabels = ["across laminae", "within laminae"];
set(gca, 'fontsize', 14);

sig_star = zeros(1,clust_data.num_clusters);
for i = 1:clust_data.num_clusters
    wb_clust_cnts = conttbl(:,i);
    wb_nclust_cnts = wb_cnts-wb_clust_cnts;
    x = table(wb_clust_cnts,wb_nclust_cnts,'VariableNames',{'inclust','ninclust'},'RowNames',{'within','between'});
    [h, p_fish, stats] = fishertest(x);
    disp("fishertest within/between cluster i: p = " + p_fish);
    sig_star(i) = p_fish < .05/clust_data.num_clusters;
end

%plot_cluster_axis(ccg_data.ccgs, clust_data.labels, clust_data.num_clusters, cmap, sig_star);
set(gcf, 'position', [314   244   563   134])
save_close_figures(flag.figure_dir + "wb_count") 

figure('position', [360   198   560   300]);
h = heatmap(round(conttbl_prop,2), 'Colormap', bone);
ax = gca();
h.XDisplayLabels = nan(size(ax.XDisplayData));
h.YDisplayLabels = ["across laminae", "within laminae"];
set(gca, 'fontsize', 14);

% construct pairwise significant test 
set(gcf, 'position', [314   244   563   134])
save_close_figures(flag.figure_dir + "wb_prop") 

%% crosstab all layer pairings, all func. ct pairings stats.
cat_var = ["cl", "sc", "ct"];
for k = 1:length(cat_var)
    pre_lab = 'pre_' + cat_var(k);
    post_lab = 'post_' + cat_var(k);
uq_layers = unique(ccg_data.(pre_lab));
uq_layers = uq_layers(~isnan(uq_layers));
pair_labels = zeros(length(clust_data.labels),1);
for i = 1:length(uq_layers)
    for j = 1:length(uq_layers)
        pair_labels(ccg_data.(pre_lab) == uq_layers(i) & ccg_data.(post_lab) == uq_layers(j)) = i + (j-1)*length(uq_layers);
    end
end

[conttbl,chi2,p,labels] = crosstab(pair_labels, clust_data.labels);
cluster_cnts = sum(conttbl,1);
wb_cnts = sum(conttbl,2);
conttbl_prop = conttbl./cluster_cnts;
disp(cat_var(k) + " x-tab all layer pairs: " + "chi^2 = " + chi2 + ", p=" + p);

for i = 1:clust_data.num_clusters
    in_clust = clust_data.labels==i;
    [conttbl_pwise,chi2_pwise,p_pwise,~] = crosstab(pair_labels, in_clust);
    disp(cat_var(k) + " x-tab all layer pairs cluster i: p = " + "chi^2 = " + chi2_pwise + ", p=" + p_pwise);
end
end

%% crosstab p-matrix, undirected
cat_var = ["cl", "sc", "ct"];
cat_lab = ["layer pairs", "func pairs", "type pairs"];
for k = 1:length(cat_var)
    pre_lab = 'pre_' + cat_var(k);
    post_lab = 'post_' + cat_var(k);
    uq_layers = unique(ccg_data.(pre_lab));
    uq_layers = uq_layers(~isnan(uq_layers));
    pair_labels = zeros(length(clust_data.labels),1);
    for i = 1:length(uq_layers)
        for j = i:length(uq_layers)
            idx = (ccg_data.(pre_lab) == uq_layers(i) & ccg_data.(post_lab) == uq_layers(j)) | (ccg_data.(post_lab) == uq_layers(i) & ccg_data.(pre_lab) == uq_layers(j));          
            pair_labels(idx) = i + (j-1)*length(uq_layers);
        end
    end
    for i = 1:clust_data.num_clusters
        for j = 1:clust_data.num_clusters
            in_clust_i = clust_data.labels==i;
            in_clust_j = clust_data.labels==j;
            shared_idx = in_clust_i | in_clust_j;
            [conttbl_pwise,chi2_pwise(i,j),p_pwise(i,j),~] = crosstab(pair_labels(shared_idx), in_clust_i(shared_idx));
            %disp(cat_var(k) + " x-tab all layer pairs cluster i: p = " + "chi^2 = " + chi2_pwise + ", p=" + p_pwise);
        end
    end
    figure('position', [360   314   232   219]);
    vals = round(log10(p_pwise));
    h = heatmap(vals,'YDisplayLabels', {'','','','',''}, 'XDisplayLabels', {'','','','',''},'MissingDataColor', 'w', 'GridVisible', 'on', 'MissingDataLabel', " ");
    colormap(bone);
    colorbar('off')
    title("log10(p), " + cat_lab(k) + " undir.");
    set(gca, 'FontName', 'Helvetica', 'FontSize', 10);
    save_close_figures(flag.figure_dir + "clust_pmtx_undir_" + cat_lab(k)) 
end

%% crosstab p-matrix, undirected within vs beteween
cat_var = ["cl", "sc", "ct"];
label_var = ["same vs diff layer", "same vs diff func. type", "same vs diff cell type"];
for k = 1:1
    pre_lab = 'pre_' + cat_var(k);
    post_lab = 'post_' + cat_var(k);
    uq_layers = unique(ccg_data.(pre_lab));
    uq_layers = uq_layers(~isnan(uq_layers));
    pair_labels = zeros(length(clust_data.labels),1);
    
    idx = (ccg_data.(pre_lab) == ccg_data.(post_lab));
    pair_labels(idx) = 1;
    
    for i = 1:clust_data.num_clusters
        for j = 1:clust_data.num_clusters
            in_clust_i = clust_data.labels==i;
            in_clust_j = clust_data.labels==j;
            shared_idx = in_clust_i | in_clust_j;
            [conttbl_pwise,chi2_pwise(i,j),p_pwise(i,j),~] = crosstab(pair_labels(shared_idx), in_clust_i(shared_idx));
            %disp(cat_var(k) + " x-tab all layer pairs cluster i: p = " + "chi^2 = " + chi2_pwise + ", p=" + p_pwise);
        end
    end
    figure('position', [360   314   232   219]);
    vals = round(log10(p_pwise),2,'significant');
    h = heatmap(vals,'YDisplayLabels', {'','','','',''}, 'XDisplayLabels', {'','','','',''},'MissingDataColor', 'w', 'GridVisible', 'on', 'MissingDataLabel', " ");
    colormap(bone);
    colorbar('off');
    set(gca, 'FontName', 'Helvetica', 'FontSize', 10);
    title("log10(p), " + label_var(k) + " undir.");
    save_close_figures(flag.figure_dir + "clust_pmtx_win_between_undir_" + label_var(k)) 
end

%% crosstab p-matrix, directed
cat_var = ["cl", "sc", "ct"];
for k = 1:length(cat_var)
    pre_lab = 'pre_' + cat_var(k);
    post_lab = 'post_' + cat_var(k);
    uq_layers = unique(ccg_data.(pre_lab));
    uq_layers = uq_layers(~isnan(uq_layers));
    pair_labels = zeros(length(clust_data.labels),1);
    for i = 1:length(uq_layers)
        for j = 1:length(uq_layers)
            idx = (ccg_data.(pre_lab) == uq_layers(i) & ccg_data.(post_lab) == uq_layers(j));          
            pair_labels(idx) = i + (j-1)*length(uq_layers);
        end
    end
    for i = 1:clust_data.num_clusters
        for j = 1:clust_data.num_clusters
            in_clust_i = clust_data.labels==i;
            in_clust_j = clust_data.labels==j;
            shared_idx = in_clust_i | in_clust_j;
            [conttbl_pwise,chi2_pwise(i,j),p_pwise(i,j),~] = crosstab(pair_labels(shared_idx), in_clust_i(shared_idx));
            %disp(cat_var(k) + " x-tab all layer pairs cluster i: p = " + "chi^2 = " + chi2_pwise + ", p=" + p_pwise);
        end
    end
    figure('position', [360   314   232   219]);
    vals = round(log10(p_pwise),2,'significant');
    h = heatmap(vals,'YDisplayLabels', {'','','','',''}, 'XDisplayLabels', {'','','','',''},'MissingDataColor', 'w', 'GridVisible', 'on', 'MissingDataLabel', " ");
    colormap(bone);
    colorbar('off');
    set(gca, 'FontName', 'Helvetica', 'FontSize', 10);
    title("log10(p), " + cat_lab(k) + " dir.");
    save_close_figures(flag.figure_dir + "clust_pmtx_dir_" + cat_var(k)) 
end

%% undirected tirin idea
cat_var = ["cl", "sc", "ct"];
poses = [   292   198   850   420;     369    72   191   420; 562    71   327   420];
for k = 1:length(cat_var)
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
                in_layer_pair = (ccg_data.(pre_lab) == uq_layers(i) & ccg_data.(post_lab) == uq_layers(j)) |(ccg_data.(pre_lab) == uq_layers(j) & ccg_data.(post_lab) == uq_layers(i));
                in_cluster = clust_data.labels == kk;
                xx = sum(in_layer_pair&in_cluster, 'all');
                pair_labels(cnt,kk) =  xx;
            end
            pair_lay_labels{cnt} = cur_labels(uq_layers(i)) + ", " + cur_labels(uq_layers(j));
        end
    end
    
    figure('position', poses(k,:));
    x = 1:.1:length(pair_lay_labels);
    vq = [];
    xs = [];
    for i = 1:clust_data.num_clusters
        vq(:,i) = interp1(1:length(pair_lay_labels), pair_labels(:,i), x, 'pchip');
        xs(:,i) = x;
    end
    p = plot(xs, vq,'linewidth', 2);
    for i = 1:clust_data.num_clusters
        p(i).Color = cmap(i,:);
    end
    set(gca, 'xticklabels', pair_lay_labels, 'xtick', 1:size(pair_labels,1));
    xtickangle(90);
    currx = xlim;
    ylabel("number of pairs");  
    xlim([0.5, currx(2)+.5])
    set_axis_defaults();
    save_close_figures(flag.figure_dir + "clust_prop_layer" + int2str(i)) 
end

%% undirected surface idea
cat_var = ["cl", "sc", "ct"];
poses = [360   195   492   423;       360   390   336   228;  360   340   321   278];

for k = 1:length(cat_var)
    pre_lab = 'pre_' + cat_var(k);
    post_lab = 'post_' + cat_var(k);
    
    uq_layers = unique(ccg_data.(pre_lab));
    uq_layers = uq_layers(~isnan(uq_layers));
    cur_labels = ccg_data.(cat_var(k) + "_labels");
    pair_labels  = [];
    pair_lay_labels = [];
    cnt = 0;
    for i = 1:length(uq_layers)
        for j = 1:length(uq_layers)
            cnt = cnt + 1;
            for kk = 1:clust_data.num_clusters
                in_layer_pair = (ccg_data.(pre_lab) == uq_layers(i) & ccg_data.(post_lab) == uq_layers(j)) | (ccg_data.(pre_lab) == uq_layers(j) & ccg_data.(post_lab) == uq_layers(i));
                in_cluster = clust_data.labels == kk;
                xx = sum(in_layer_pair&in_cluster, 'all');
                pair_labels(i,j,kk) =  xx;
            end
            pair_lay_labels{i,j} = cur_labels(uq_layers(i)) + ", " + cur_labels(uq_layers(j));
        end
    end
    for vs = 1:2
    figure('position', poses(k,:));
    [x,y] = meshgrid(1:.1:length(pair_lay_labels));
    [xo,yo] = meshgrid(1:length(uq_layers));
    vq = [];
    xs = [];
    cc = [];
    for i = 1:clust_data.num_clusters
        vq= interp2(xo,yo, pair_labels(:,:,i), x,y, 'cubic');
        %xs(:,i) = x;
        for aa = 1:size(x*y)
            for bb = 1:size(x*y)
                cc(aa,bb,:) = cmap(i,:);
            end
        end
        s = surf(x,y,vq, cc, 'linestyle', 'none');
        hold on;
    end

    %p = plot(xs, vq,'linewidth', 2);
    %for i = 1:clust_data.num_clusters
   %    p(i).Color = cmap(i,:);
    %end
    %surf(xs,xs,vq
    set(gca, 'xticklabels', cur_labels, 'xtick', 1:size(pair_labels,1),  'yticklabels', cur_labels, 'ytick', 1:size(pair_labels,1)); 
    zlabel("number of pairs");
    %xtickangle(90);
    currx = xlim;
    %xlim([0.5, currx(2)+.5])
    set_axis_defaults();
    set(gca, 'fontsize', 12);
    grid off;
    if k<3 && vs == 1
        set(gca, 'view',[  135   40]);
    elseif vs == 2
        view(2);
    end
    save_close_figures(flag.figure_dir + "undir_clust_surf_layer_" + cat_var(k) + int2str(vs))
    end
end

%% directed surface idea
cat_var = ["cl", "sc", "ct"];
poses = [360   195   492   423;       360   390   336   228;  360   340   321   278];

for k = 1:length(cat_var)
    pre_lab = 'pre_' + cat_var(k);
    post_lab = 'post_' + cat_var(k);
    
    uq_layers = unique(ccg_data.(pre_lab));
    uq_layers = uq_layers(~isnan(uq_layers));
    cur_labels = ccg_data.(cat_var(k) + "_labels");
    pair_labels  = [];
    pair_lay_labels = [];
    cnt = 0;
    for i = 1:length(uq_layers)
        for j = 1:length(uq_layers)
            cnt = cnt + 1;
            for kk = 1:clust_data.num_clusters
                in_layer_pair = (ccg_data.(pre_lab) == uq_layers(i) & ccg_data.(post_lab) == uq_layers(j));
                in_cluster = clust_data.labels == kk;
                xx = sum(in_layer_pair&in_cluster, 'all');
                pair_labels(i,j,kk) =  xx;
            end
            pair_lay_labels{i,j} = cur_labels(uq_layers(i)) + ", " + cur_labels(uq_layers(j));
        end
    end
    for vs = 1:2
    figure('position', poses(k,:));
    [x,y] = meshgrid(1:.1:length(pair_lay_labels));
    [xo,yo] = meshgrid(1:length(uq_layers));
    vq = [];
    xs = [];
    cc = [];
    for i = 1:clust_data.num_clusters
        vq= interp2(xo,yo, pair_labels(:,:,i), x,y, 'cubic');
        %xs(:,i) = x;
        for aa = 1:size(x*y)
            for bb = 1:size(x*y)
                cc(aa,bb,:) = cmap(i,:);
            end
        end
        s = surf(x,y,vq, cc, 'linestyle', 'none');
        hold on;
    end

    %p = plot(xs, vq,'linewidth', 2);
    %for i = 1:clust_data.num_clusters
   %    p(i).Color = cmap(i,:);
    %end
    %surf(xs,xs,vq
    set(gca, 'xticklabels', cur_labels, 'xtick', 1:size(pair_labels,1),  'yticklabels', cur_labels, 'ytick', 1:size(pair_labels,1)); 
    xlabel("lag");
    ylabel("lead");
    zlabel("number of pairs");
    %xtickangle(90);
    currx = xlim;
    %xlim([0.5, currx(2)+.5])
    set_axis_defaults();
    set(gca, 'fontsize', 12);
    grid off;
    if k<3 && vs == 1
        set(gca, 'view',[  135   40]);
    elseif vs == 2
        view(2);
    end
    save_close_figures(flag.figure_dir + "clust_surf_layer_" + cat_var(k) + int2str(vs))
    end
end

%% directed tirin idea
cat_var = ["cl", "sc", "ct"];
for k = 1:length(cat_var)
    pre_lab = 'pre_' + cat_var(k);
    post_lab = 'post_' + cat_var(k);
uq_layers = unique(ccg_data.(pre_lab));
uq_layers = uq_layers(~isnan(uq_layers));
pair_labels = zeros(length(clust_data.labels),1);
for i = 1:length(uq_layers)
    for j = 1:length(uq_layers)
        pair_labels(ccg_data.(pre_lab) == uq_layers(i) & ccg_data.(post_lab) == uq_layers(j)) = i + (j-1)*length(uq_layers);
    end
end

[conttbl,chi2,p,labels] = crosstab(pair_labels, clust_data.labels);
cluster_cnts = sum(conttbl,1);
wb_cnts = sum(conttbl,2);
conttbl_prop = conttbl./cluster_cnts;
disp(cat_var(k) + " x-tab all layer pairs: " + "chi^2 = " + chi2 + ", p=" + p);

for i = 1:clust_data.num_clusters
    in_clust = clust_data.labels==i;
    [conttbl_pwise,chi2_pwise,p_pwise,~] = crosstab(pair_labels, in_clust);
    disp(cat_var(k) + " x-tab all layer pairs cluster i: p = " + "chi^2 = " + chi2_pwise + ", p=" + p_pwise);
end
end

%% undirected within/between tirin idea
cat_var = ["cl", "sc", "ct"];
cat_labs = ["laminae", "func type", "type"];
poses = [562    71   327   420;    562    71   327   420; 562    71   327   420];
for k = 1:length(cat_var)
    pre_lab = 'pre_' + cat_var(k);
    post_lab = 'post_' + cat_var(k);
    
    uq_layers = unique(ccg_data.(pre_lab));
    uq_layers = uq_layers(~isnan(uq_layers));
    cur_labels = ccg_data.(cat_var(k) + "_labels");
    pair_labels  = [];
    pair_lay_labels = [];
    cnt = 0;
    
    cnt = cnt + 1;
    for kk = 1:clust_data.num_clusters
        in_layer_pair = (ccg_data.(pre_lab) == ccg_data.(post_lab));
        in_cluster = clust_data.labels == kk;
        xx = sum(in_layer_pair&in_cluster, 'all');
        pair_labels(cnt,kk) =  xx;
    end
    pair_lay_labels{cnt} = "within " + cat_labs(k);
 
    cnt = cnt + 1;
    for kk = 1:clust_data.num_clusters
        in_layer_pair = (ccg_data.(pre_lab) ~= ccg_data.(post_lab));
        in_cluster = clust_data.labels == kk;
        xx = sum(in_layer_pair&in_cluster, 'all');
        pair_labels(cnt,kk) =  xx;
    end
    pair_lay_labels{cnt} = "between " + cat_labs(k);
    
    figure('position', poses(k,:));
    x = 1:.1:length(pair_lay_labels);
    vq = [];
    xs = [];
    for i = 1:clust_data.num_clusters
        vq(:,i) = interp1(1:length(pair_lay_labels), pair_labels(:,i), x, 'pchip');
        xs(:,i) = x;
    end
    p = plot(xs, vq,'linewidth', 2);
    for i = 1:clust_data.num_clusters
        p(i).Color = cmap(i,:);
    end
    set(gca, 'xticklabels', pair_lay_labels, 'xtick', 1:size(pair_labels,1));
    xtickangle(90);
    currx = xlim;
    ylabel("number of pairs");  
    xlim([0.5, currx(2)+.5])
    set_axis_defaults();
    save_close_figures(flag.figure_dir + "clust_prop_withinbet" + int2str(k)) 
end

%% undirected bar plots
cat_var = ["cl", "sc", "ct"];
cat_labs = ["laminae", "func type", "type"];
for jj = 1:length(cat_var)
    pre_lab = "pre_" + cat_var(jj);
    post_lab = "post_" + cat_var(jj);
    labs = cat_var(jj) + "_labels";
    for i = 1:length(unique(clust_data.labels))
        cluster_size(i) = sum(clust_data.labels == i);
        for pre_cl = 1:length(ccg_data.(labs))
            for post_cl = 1:length(ccg_data.(labs))
                cnt_cler_pair(i, pre_cl, post_cl) = sum((ccg_data.(pre_lab)(clust_data.labels ==i) == pre_cl & ccg_data.(post_lab)(clust_data.labels ==i) == post_cl) |...
                    (ccg_data.(pre_lab)(clust_data.labels ==i) == post_cl & ccg_data.(post_lab)(clust_data.labels ==i) == pre_cl));
                norm_cler_pair(i, pre_cl, post_cl) = cnt_cler_pair(i, pre_cl, post_cl)/cluster_size(i);
           end
        end
    end
    for pre_cl = 1:length(ccg_data.(labs))
        for post_cl = 1:length(ccg_data.(labs))
            if post_cl<=pre_cl
            figure('position', [  360.0000  445.0000  240.3333  173.0000]);
            b = bar(1:clust_data.num_clusters,norm_cler_pair(:,pre_cl, post_cl),'facecolor', 'flat', 'CData', cmap, 'edgecolor', 'w', 'linewidth', 2);

            hold on;
            %errorbar(1:clust_data.num_clusters, at_mean.(attributes{k}), at_sem.(attributes{k}),'.k', 'linewidth', 1.2, 'markersize', 1);
            set_axis_defaults();
            set(gca, 'fontsize', 10);
            ylab = ccg_data.(labs)(pre_cl) + " + " + ccg_data.(labs)(post_cl) + " cnt/clust. sz";
            ylabel(ylab, 'fontsize', 10);
            set(gca, 'xticklabel', {}, 'xtick', []);
            xlim([.25, clust_data.num_clusters + .75])
            set_axis_defaults();

            figure('position', [  360.0000  445.0000  240.3333  173.0000]);
            b = bar(1:clust_data.num_clusters,cnt_cler_pair(:,pre_cl, post_cl),'facecolor', 'flat', 'CData', cmap, 'edgecolor', 'w', 'linewidth', 2);

            hold on;
            %errorbar(1:clust_data.num_clusters, at_mean.(attributes{k}), at_sem.(attributes{k}),'.k', 'linewidth', 1.2, 'markersize', 1);
            set_axis_defaults();
            set(gca, 'fontsize', 10);
            ylabel(ccg_data.(labs)(pre_cl) + " + " + ccg_data.(labs)(post_cl) + " cnt" , 'fontsize', 10);
            set(gca, 'xticklabel', {}, 'xtick', []);
            xlim([.25, clust_data.num_clusters + .75])
            set_axis_defaults();

            %overall chi2 test
            in_pair = (ccg_data.(pre_lab) == pre_cl & ccg_data.(post_lab) == post_cl)|(ccg_data.(pre_lab) == post_cl & ccg_data.(post_lab) == pre_cl);
            [tbl, chi2, p] = crosstab(clust_data.labels, in_pair);
            if p<.05/60
                text(1,1,'*', 'units', 'normalized', 'fontsize', 20);
            end

            %pairwise chi2 test
            y = cnt_cler_pair(:,pre_cl, post_cl);
            [m, i] = max(y);
            for c_cnt = 1:clust_data.num_clusters
                 sub_idx = clust_data.labels == i | clust_data.labels == c_cnt;
                 [tbl, chi2, p] = crosstab(clust_data.labels(sub_idx), in_pair(sub_idx));
                 if p<.05/60
                    text(c_cnt,cnt_cler_pair(c_cnt,pre_cl, post_cl)+2,'*', 'units', 'data','horizontalalignment', 'center','fontsize', 14);
                 end
            end
            save_close_figures(flag.figure_dir +erase(erase(erase(ylab, '/'),'.'),'\'));
            end
        end
    end
end

%% directed bar plots
cat_var = ["cl", "sc", "ct"];
cat_labs = ["laminae", "func type", "type"];
for jj = 1:length(cat_var)
    pre_lab = "pre_" + cat_var(jj);
    post_lab = "post_" + cat_var(jj);
    labs = cat_var(jj) + "_labels";
    for i = 1:length(unique(clust_data.labels))
        cluster_size(i) = sum(clust_data.labels == i);
        for pre_cl = 1:length(ccg_data.(labs))
            for post_cl = 1:length(ccg_data.(labs))
                cnt_cler_pair(i, pre_cl, post_cl) = sum((ccg_data.(pre_lab)(clust_data.labels ==i) == pre_cl & ccg_data.(post_lab)(clust_data.labels ==i) == post_cl));
                norm_cler_pair(i, pre_cl, post_cl) = cnt_cler_pair(i, pre_cl, post_cl)/cluster_size(i);
            end
        end
    end
    for pre_cl = 1:length(ccg_data.(labs))
        for post_cl = 1:length(ccg_data.(labs))
            figure('position', [  360.0000  445.0000  240.3333  173.0000]);
            b = bar(1:clust_data.num_clusters,norm_cler_pair(:,pre_cl, post_cl),'facecolor', 'flat', 'CData', cmap, 'edgecolor', 'w', 'linewidth', 2);

            hold on;
            %errorbar(1:clust_data.num_clusters, at_mean.(attributes{k}), at_sem.(attributes{k}),'.k', 'linewidth', 1.2, 'markersize', 1);
            set_axis_defaults();
            set(gca, 'fontsize', 10);
            ylab = ccg_data.(labs)(pre_cl) + " to " + ccg_data.(labs)(post_cl) + " cnt/clust. sz";
            ylabel(ylab, 'fontsize', 10);
            set(gca, 'xticklabel', {}, 'xtick', []);
            xlim([.25, clust_data.num_clusters + .75])
            set_axis_defaults();

            figure('position', [360.0000  445.0000  240.3333  173.0000]);
            b = bar(1:clust_data.num_clusters,cnt_cler_pair(:,pre_cl, post_cl),'facecolor', 'flat', 'CData', cmap, 'edgecolor', 'w', 'linewidth', 2);

            hold on;
            %errorbar(1:clust_data.num_clusters, at_mean.(attributes{k}), at_sem.(attributes{k}),'.k', 'linewidth', 1.2, 'markersize', 1);
            set_axis_defaults();
            set(gca, 'fontsize', 10);
            ylabel(ccg_data.(labs)(pre_cl) + " to " + ccg_data.(labs)(post_cl) + " cnt" , 'fontsize', 10);
            set(gca, 'xticklabel', {}, 'xtick', []);
            xlim([.25, clust_data.num_clusters + .75])
            set_axis_defaults();
            
            %overall chi2 test
            in_pair = (ccg_data.(pre_lab) == pre_cl & ccg_data.(post_lab) == post_cl);
            [tbl, chi2, p] = crosstab(clust_data.labels, in_pair);
            if p<.05/60
                text(1,1,'*', 'units', 'normalized', 'fontsize', 20);
            end

            %pairwise chi2 test
            y = cnt_cler_pair(:,pre_cl, post_cl);
            [m, i] = max(y);
            for c_cnt = 1:clust_data.num_clusters
                 sub_idx = clust_data.labels == i | clust_data.labels == c_cnt;
                 [tbl, chi2, p] = crosstab(clust_data.labels(sub_idx), in_pair(sub_idx));
                 if p<.05/60
                    text(c_cnt,cnt_cler_pair(c_cnt,pre_cl, post_cl)+2,'*', 'units', 'data','horizontalalignment', 'center','fontsize', 14);
                 end
            end
            save_close_figures(flag.figure_dir +"zz" + erase(erase(erase(ylab, '/'),'.'),'\'));
        end
    end
end


%% directed bar plot intersections
cat_var = ["sc_x_cl"];
cat_labs = ["laminae x func type"];
pre_cl_is_4 = ccg_data.pre_cl == 2 | ccg_data.pre_cl == 3 | ccg_data.pre_cl == 4 ;
post_cl_is_23 = ccg_data.post_cl == 1;
pre_cl_is_23 = ccg_data.pre_cl == 1 ;
post_cl_is_4 = ccg_data.post_cl == 2 | ccg_data.post_cl == 3 | ccg_data.post_cl == 4 ;

lay_labs = (pre_cl_is_4 + pre_cl_is_23*2);
tot_labs = lay_labs+3*(ccg_data.pre_sc-1);
ccg_data.pre_sc_x_cl = tot_labs+1;

lay_labs = (post_cl_is_4 + post_cl_is_23*2);
tot_labs = lay_labs+3*(ccg_data.post_sc-1);
ccg_data.post_sc_x_cl = tot_labs+1;

ccg_data.sc_x_cl_labels = ["not4 or 23,comp","4,comp", "2/3,comp","not4 or 23,simp","4,simp", "2/3,simp"];

jj = 1;
pre_lab = "pre_" + cat_var(jj);
post_lab = "post_" + cat_var(jj);
labs = cat_var(jj) + "_labels";
for i = 1:length(unique(clust_data.labels))
    cluster_size(i) = sum(clust_data.labels == i);
    for pre_cl = 1:length(ccg_data.(labs))
        for post_cl = 1:length(ccg_data.(labs))
            cnt_cler_pair(i, pre_cl, post_cl) = sum((ccg_data.(pre_lab)(clust_data.labels ==i) == pre_cl & ccg_data.(post_lab)(clust_data.labels ==i) == post_cl));
            norm_cler_pair(i, pre_cl, post_cl) = cnt_cler_pair(i, pre_cl, post_cl)/cluster_size(i);
        end
    end
end
for pre_cl = 1:length(ccg_data.(labs))
    for post_cl = 1:length(ccg_data.(labs))
        if pre_cl ~= 1 && pre_cl ~= 4 && post_cl~=1 && post_cl ~=4
            figure('position', [  360.0000  445.0000  240.3333  173.0000]);
            b = bar(1:clust_data.num_clusters,norm_cler_pair(:,pre_cl, post_cl),'facecolor', 'flat', 'CData', cmap, 'edgecolor', 'w', 'linewidth', 2);
            
            hold on;
            %errorbar(1:clust_data.num_clusters, at_mean.(attributes{k}), at_sem.(attributes{k}),'.k', 'linewidth', 1.2, 'markersize', 1);
            set_axis_defaults();
            set(gca, 'fontsize', 10);
            ylab = ccg_data.(labs)(pre_cl) + " to " + ccg_data.(labs)(post_cl) + " cnt/clust. sz";
            ylabel(ylab, 'fontsize', 10);
            set(gca, 'xticklabel', {}, 'xtick', []);
            xlim([.25, clust_data.num_clusters + .75])
            set_axis_defaults();
            
            figure('position', [360.0000  445.0000  240.3333  173.0000]);
            b = bar(1:clust_data.num_clusters,cnt_cler_pair(:,pre_cl, post_cl),'facecolor', 'flat', 'CData', cmap, 'edgecolor', 'w', 'linewidth', 2);
            
            hold on;
            %errorbar(1:clust_data.num_clusters, at_mean.(attributes{k}), at_sem.(attributes{k}),'.k', 'linewidth', 1.2, 'markersize', 1);
            set_axis_defaults();
            set_axis_defaults();
            set(gca, 'fontsize', 8);
            ylabel(ccg_data.(labs)(pre_cl) + " to " + ccg_data.(labs)(post_cl) + " cnt" , 'fontsize', 10);
            set(gca, 'xticklabel', {}, 'xtick', []);
            xlim([.25, clust_data.num_clusters + .75])
            
            %overall chi2 test
            in_pair = (ccg_data.(pre_lab) == pre_cl & ccg_data.(post_lab) == post_cl);
            [tbl, chi2, p] = crosstab(clust_data.labels, in_pair);
            if p<.05/60
                text(1,1,'*', 'units', 'normalized', 'fontsize', 20);
            end
            
            %pairwise chi2 test
            y = cnt_cler_pair(:,pre_cl, post_cl);
            [m, i] = max(y);
            for c_cnt = 1:clust_data.num_clusters
                sub_idx = clust_data.labels == i | clust_data.labels == c_cnt;
                [tbl, chi2, p] = crosstab(clust_data.labels(sub_idx), in_pair(sub_idx));
                if p<.05/60
                    text(c_cnt,cnt_cler_pair(c_cnt,pre_cl, post_cl)+2,'*', 'units', 'data','horizontalalignment', 'center','fontsize', 14);
                end
            end
            save_close_figures(flag.figure_dir +"zz" + erase(erase(erase(ylab, '/'),'.'),'\'));
        end
    end
end


%% time to fs
figure('position', [  360.0000  445.0000  240.3333  173.0000]);
for i = 1:clust_data.num_clusters
    time_to_fs = (ccg_data.cluster.time_to_fs(ccg_data.pre_id(clust_data.labels == i)) +ccg_data.cluster.time_to_fs(ccg_data.post_id(clust_data.labels==i)));
    mtime_to_fs(i) = nanmean(time_to_fs);
    semtime_to_fs(i) = nansem(time_to_fs);
    hold on;
end

b = bar(1:clust_data.num_clusters,mtime_to_fs,'facecolor', 'flat', 'CData', cmap, 'edgecolor', 'w', 'linewidth', 2);

hold on;
errorbar(1:clust_data.num_clusters, mtime_to_fs, semtime_to_fs,'.k', 'linewidth', 1.2, 'markersize', 1);
set_axis_defaults();
set(gca, 'fontsize', 10);
ylab = "t. to first spike sum";
ylabel(ylab , 'fontsize', 10);
set(gca, 'xticklabel', {}, 'xtick', []);
xlim([.25, clust_data.num_clusters + .75])
ylim([min(mtime_to_fs)-.1*min(mtime_to_fs), max(mtime_to_fs)+.1*min(mtime_to_fs)]);
set_axis_defaults();
save_close_figures(flag.figure_dir +erase(erase(erase(ylab, '/'),'.'),'\'));


figure('position', [  360.0000  445.0000  240.3333  173.0000]);
for i = 1:clust_data.num_clusters
    time_to_fs = abs(ccg_data.cluster.time_to_fs(ccg_data.pre_id(clust_data.labels == i))-ccg_data.cluster.time_to_fs(ccg_data.post_id(clust_data.labels==i)));
    mtime_to_fs(i) = nanmean(time_to_fs);
    semtime_to_fs(i) = nansem(time_to_fs);
    hold on;
end

b = bar(1:clust_data.num_clusters,mtime_to_fs,'facecolor', 'flat', 'CData', cmap, 'edgecolor', 'w', 'linewidth', 2);

hold on;
errorbar(1:clust_data.num_clusters, mtime_to_fs, semtime_to_fs,'.k', 'linewidth', 1.2, 'markersize', 1);
set_axis_defaults();
set(gca, 'fontsize', 10);
ylab = "t. to first spike diff";
ylabel(ylab, 'fontsize', 10);
set(gca, 'xticklabel', {}, 'xtick', []);
xlim([.25, clust_data.num_clusters + .75])
ylim([min(mtime_to_fs)-.1*min(mtime_to_fs), max(mtime_to_fs)+.1*min(mtime_to_fs)]);
set_axis_defaults();
save_close_figures(flag.figure_dir +erase(erase(erase(ylab, '/'),'.'),'\'));

%% heatmaps prop across clusters
figure('position', [  377.0000  198.3333  511.3333  408.0000]);
tbl = crosstab(ccg_data.pre_cl(clust_data.labels ==i), ccg_data.post_cl(clust_data.labels ==i));
incolor = cmap(i,:);
scale = linspace(1,0,100);
shades = double(incolor.*scale.');
tints = incolor.*double(uint8(incolor + (255-incolor).*fliplr(scale).'))./255;
cont_heatmap(tbl, ccg_data.cl_labels, colormap(tints));
save_close_figures(flag.figure_dir + "clust_cnt_layer" + int2str(i))
for i = 1:length(unique(clust_data.labels))
    figure('position', [  377.0000  198.3333  511.3333  408.0000]);
    tbl = crosstab(ccg_data.pre_sc(clust_data.labels ==i), ccg_data.post_sc(clust_data.labels ==i));
    incolor = cmap(i,:);
    scale = linspace(1,0,100);
    shades = double(incolor.*scale.');
    tints = incolor.*double(uint8(incolor + (255-incolor).*fliplr(scale).'))./255;
    cont_heatmap(tbl, ccg_data.sc_labels, colormap(tints));
    save_close_figures(flag.figure_dir + "clust_cnt_sc" + int2str(i))
end