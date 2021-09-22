function compute_clusters()
addpath(genpath(pwd));
flag = config();

for group_idx = 1:flag.group_cnt
% load postprocessing output
load(flag.postproc_output(group_idx), 'ccg_data');

% subset ccgs based on significance criteria
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
    
disp("total number of pairs: " + length(sig_idx));
disp("number of significant pairs: " + sum(sig_idx));

% setup clustering
rng(flag.rng_seed);

X = ccg_data.ccgs;

% dimen. reduce
Y = tsne(X, 'NumDimensions', flag.tsne_dims); 

% kmeans
clusters = zeros(size(Y,1),flag.max_k);
wss = zeros(flag.max_k, 1);
for i = 1:flag.max_k
    clusters(:,i) = kmeans(Y,i,'replicates',flag.kmeans_rep, 'maxiter', flag.kmeans_maxiter);
    wss(i) = within_cluster_ss(Y, clusters(:,i));
end

sil_eva = evalclusters(Y,clusters,'silhouette');

% setup output structure
clust_data.labels = clusters(:,sil_eva.OptimalK);
clust_data.num_clusters = sil_eva.OptimalK;
clust_data.clusters = clusters;
clust_data.sil_eva = sil_eva;
clust_data.tsne_mtx = Y;
clust_data.wss = wss;
clust_data.clust_flag = flag;

% save output structure
save(flag.cluster_output(group_idx), 'clust_data');
end