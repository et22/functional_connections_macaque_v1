function [ses1, ses2, ses3, layer_centers] = scale_depth_by_layer(ses3, ses2, ses1)
layers = unique(ses3.Cluster_celllayer);
layers = layers(~isnan(layers));
for i = 1:length(layers)
    layer_boundaries(i,1)  = min(ses3.Cluster_celldepth(layers(i) == ses3.Cluster_celllayer));
    layer_boundaries(i,2)  = max(ses3.Cluster_celldepth(layers(i) == ses3.Cluster_celllayer));    
    layer_centers(i) = mean(layer_boundaries(i,:));
    ses2.Cluster_celldepth(layers(i) == ses2.Cluster_celllayer) = rescale(ses2.Cluster_celldepth(layers(i) == ses2.Cluster_celllayer), layer_boundaries(i,1), layer_boundaries(i,2));
    ses1.Cluster_celldepth(layers(i) == ses1.Cluster_celllayer) = rescale(ses1.Cluster_celldepth(layers(i) == ses1.Cluster_celllayer), layer_boundaries(i,1), layer_boundaries(i,2));
end

ses3.Cluster_celldepth_scale = ses3.Cluster_celldepth;
ses2.Cluster_celldepth_scale = ses3.Cluster_celldepth;
ses1.Cluster_celldepth_scale = ses3.Cluster_celldepth;

end