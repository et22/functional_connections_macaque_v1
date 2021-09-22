function idx = sort_clusters_by_lag(num_clusters, smoothed_ccgs, idx)
for i = 1:num_clusters
    [~,lags] = max(smoothed_ccgs(idx==i,:),[],2);
    lag(i) = mean(abs(lags-11));
end
[b, idxes] = sort(lag);
temp_idx = idx;
for i = 1:num_clusters
    idx(temp_idx == idxes(i)) = i; 
end