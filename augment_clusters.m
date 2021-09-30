%function augment_clusters()
addpath(genpath(pwd));
flag = config();

% load augment ccgs and subset based on significance criteria
load(flag.augment_augment_path, 'ccg_data');
ccg_data = ccg_data.ccg;
[aug_ccg_data, sig_idx] = get_significant_ccgs(ccg_data, flag);
disp("total number of pairs: " + length(sig_idx));
disp("number of significant pairs: " + sum(sig_idx));

% load base ccgs and subset based on significance criteria
load(flag.augment_base_path, 'ccg_data');
ccg_data = ccg_data.ccg;
[base_ccg_data, sig_idx] = get_significant_ccgs(ccg_data, flag);
disp("total number of pairs: " + length(sig_idx));
disp("number of significant pairs: " + sum(sig_idx));

% load clustered ccgs to augment
load(flag.augment_cluster_path, 'clust_data');

augmented_ccgs = aug_ccg_data.ccgs;
clustered_ccgs = base_ccg_data.ccgs;

augment_labels = zeros(length(augmented_ccgs), 1);
for augment_ccg = 1:length(augmented_ccgs)
    ccg_aug = augmented_ccgs(augment_ccg,:);
    dist = zeros(1,length(clustered_ccgs));
    for cluster_ccg = 1:length(clustered_ccgs)
        ccg_clust = clustered_ccgs(cluster_ccg,:);
        dist(cluster_ccg) = sqrt(sum((ccg_aug - ccg_clust) .^ 2));
    end
    clust_dist = zeros(1,clust_data.num_clusters);
    
    for clust_idx = 1:clust_data.num_clusters
        clust_dist(clust_idx) = median(dist(clust_data.labels == clust_idx));
    end
    [~,augment_labels(augment_ccg)] = min(clust_dist);
end

clust_data.labels = [clust_data.labels; augment_labels];

save(flag.augment_output, 'clust_data');
%end

