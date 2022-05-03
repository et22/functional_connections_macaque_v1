function compute_clusters()
addpath(genpath(pwd));
flag = config();

for group_idx = 1:flag.group_cnt
% load postprocessing output
load(flag.postproc_output(group_idx), 'ccg_data');

% subset ccgs based on significance criteria
ccg_data = ccg_data.ccg;
[ccg_data, sig_idx] = get_significant_ccgs(ccg_data, flag);
disp("total number of pairs: " + length(sig_idx));
disp("number of significant pairs: " + sum(sig_idx));

% setup clustering
rng(flag.rng_seed);

X = ccg_data.ccgs;

switch flag.process_cluster(flag.process_flag)
    case 'mean_center'
       X = X - nanmean(X,2); % subtract mean of each row 
       disp('mean-centering ccgs before clustering');
    case 'zscore'
       X = normalize(X,2); % zscore each row
       disp('z-scoring ccgs before clustering');
    case 'l2norm'
       X = normalize(X,'norm', 2);
       disp('l2normalizing ccgs before clustering');
end

% dimen. reduce
disp('starting tsne...');
Y = tsne(X, 'NumDimensions', flag.tsne_dims); 


% kmeans
disp('starting clustering...');
clusters = zeros(size(Y,1),flag.max_k);
wss = zeros(flag.max_k, 1);
for i = 1:flag.max_k
    clusters(:,i) = kmeans(Y,i,'replicates',flag.kmeans_rep, 'maxiter', flag.kmeans_maxiter);
    %opt.MaxIter = flag.kmeans_maxiter;
    %clusters(:,i) = kmedoids(Y,i,'replicates',flag.kmeans_rep,'Options', opt);
    
    wss(i) = within_cluster_ss(Y, clusters(:,i));
end

sil_eva = evalclusters(Y,clusters,'silhouette');

% setup output structure
if flag.fix_k
    opt_k = flag.k;
else
    opt_k = sil_eva.OptimalK;
end
clust_data.labels = clusters(:,opt_k);%sil_eva.OptimalK);
clust_data.num_clusters = opt_k;%sil_eva.OptimalK;
clust_data.clusters = clusters;
clust_data.sil_eva = sil_eva;
clust_data.tsne_mtx = Y;
clust_data.wss = wss;
clust_data.clust_flag = flag;

% save output structure
save(flag.cluster_output(group_idx), 'clust_data');
end