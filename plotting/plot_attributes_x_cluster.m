function plot_attributes_x_cluster(ccg_data,clust_data, cmap)
    attributes = ["peaks", "troughs", "peak_lag", "","","","pair_distance", "peak_width", "trough_width"];
    attributes{4} = 'r_ori';
    attributes{5} = 'r_freq';
    attributes{6} = 'r_eye';
    attribute_labels = ["peaks", "troughs", "\tau", "r_{ori}", "r_{freq}", "r_{eye}", "pair distance", "peak width", "trough width"];
    
    for i = 1:(clust_data.num_clusters)
        for k = 1:length(attributes)
            at_mean.(attributes{k})(i) = nanmean(ccg_data.(attributes{k})(clust_data.labels==i));
            at_sem.(attributes{k})(i) = nanstd(ccg_data.(attributes{k})(clust_data.labels==i))/sqrt(length(ccg_data.(attributes{k})(clust_data.labels==i)));
        end
    end
    
    for k = 1:length(attributes)
        figure('position', [  360.0000  445.0000  240.3333  173.0000]);
        b = bar(1:clust_data.num_clusters,at_mean.(attributes{k}),'facecolor', 'flat', 'CData', cmap, 'edgecolor', 'w', 'linewidth', 2);

        hold on;
        errorbar(1:clust_data.num_clusters, at_mean.(attributes{k}), at_sem.(attributes{k}),'.k', 'linewidth', 1.2, 'markersize', 1);
        set_axis_defaults();
        set(gca, 'fontsize', 10);
        ylabel(attribute_labels(k), 'fontsize', 10);
        set(gca, 'xticklabel', {}, 'xtick', []);
        xlim([.25, clust_data.num_clusters + .75])
        set_axis_defaults();
    end
end