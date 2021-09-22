function [data_included] = split_layers_at_mean(data_included, flag)
    if flag.split_layers
        depths_56 = data_included.Cluster_celldepth(data_included.Cluster_celllayer == 1);
        max_depth_56 = max(depths_56);
        min_depth_56 = min(depths_56);
        split_56 = mean([max_depth_56, min_depth_56]);
        % add space for new layers
        data_included.Cluster_celllayer = data_included.Cluster_celllayer+1; % layers shift right 
        data_included.Cluster_celllayer(data_included.Cluster_celldepth<split_56& data_included.Cluster_celllayer == 2) = 1; % lay 6
        
        depths_4cab = data_included.Cluster_celldepth(data_included.Cluster_celllayer == 3);
        max_depth_4cab = max(depths_4cab);
        min_depth_4cab = min(depths_4cab);
        split_4cab = mean([max_depth_4cab, min_depth_4cab]);
        % add space for new layers
        data_included.Cluster_celllayer( data_included.Cluster_celllayer>2) = data_included.Cluster_celllayer( data_included.Cluster_celllayer>2)+1; %>=3 -> n + 1 
        data_included.Cluster_celllayer(data_included.Cluster_celldepth<split_4cab& data_included.Cluster_celllayer == 4) = 3; % lay 4cb
    end
    
    % put WM before 6
    temp = data_included.Cluster_celllayer;
    wm_idx = max(temp);
    temp(temp==wm_idx) = 0;
    data_included.Cluster_celllayer = temp + 1;
    
    % flip array s.t. 6 is on bottom 
    layers = unique(data_included.Cluster_celllayer);
    layers = layers(~isnan(layers));
    for i = 1:floor(length(layers)/2)
        temp = data_included.Cluster_celllayer == i; 
        data_included.Cluster_celllayer(data_included.Cluster_celllayer == length(layers)-i+1) = i;
        data_included.Cluster_celllayer(temp) = length(layers)-i+1;
    end
    