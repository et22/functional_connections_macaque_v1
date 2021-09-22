function sumval = within_cluster_ss(Y, idx)
    clusters = unique(idx); 
    sumval = 0;
    for i = 1:length(clusters)
        data = Y(clusters(i)==idx,:);
        avg_clust = mean(data, 1);
        for j = 1:size(data,1)
            sumval = sumval + sum((data(j,:)-avg_clust).^2);
        end
    end
end