function [data_included] = split_layers_at_mean(data_included, flag, session_number)
ses13 =[2300 2355 2262; % max 2/3 depth
        1609 1748 1610; % max 4a/b depth
        1353 1432 1302; % max 4c depth
        1070 1190 990; % max 5/6 depth
        508 702 353; % max WM depth
        0   0   0]; % min WM depth

ses13_split =[2300 2355 2262; % max 2/3 depth
        1609 1748 1610; % max 4a/b depth
        1353 1432 1302; % max 4ca depth
        1212 1311 1146; % max 4cb depth
        1070 1190 990; % max 5 depth
        789 946  672; % max 6 depth
        508 702 353; % max WM depth
        0   0   0]; % min WM depth
    
ses45 = [0 0; % min 2/3 depth
        73 0; % max 2/3 depth
        422 312; % max 4a/b depth
        750 550; % max 4c depth
        1124 933; % max 5/6 depth
        1500 1500]; % max WM depth

ses45_split = [0 0; % min 2/3 depth
    73 0; % max 2/3 depth
    422 312; % max 4a/b depth
    586 431; % max 4ca depth
    750 550; % max 4cb depth
    937  742; % max 5 depth
    1124 933; % max 6 depth
    1500 1500]; % max WM depth

    
if session_number >= 1 && session_number <= 3
    bin_edges = ses13_split(:,session_number);
    y = discretize(-data_included.Cluster_celldepth, -bin_edges);
    data_included.Cluster_celllayer = y;
elseif session_number >= 4 && session_number <= 5
    bin_edges = ses45_split(:, session_number-3);
    y = discretize(data_included.Cluster_celldepth, bin_edges);
    data_included.Cluster_celllayer = y;
end


    % old 
%     if flag.split_layers
%         depths_56 = data_included.Cluster_celldepth(data_included.Cluster_celllayer == 1);
%         max_depth_56 = max(depths_56);
%         min_depth_56 = min(depths_56);
%         split_56 = mean([max_depth_56, min_depth_56]);
%         % add space for new layers
%         data_included.Cluster_celllayer = data_included.Cluster_celllayer+1; % layers shift right 
%         data_included.Cluster_celllayer(data_included.Cluster_celldepth<split_56& data_included.Cluster_celllayer == 2) = 1; % lay 6
%         
%         depths_4cab = data_included.Cluster_celldepth(data_included.Cluster_celllayer == 3);
%         max_depth_4cab = max(depths_4cab);
%         min_depth_4cab = min(depths_4cab);
%         split_4cab = mean([max_depth_4cab, min_depth_4cab]);
%         % add space for new layers
%         data_included.Cluster_celllayer( data_included.Cluster_celllayer>2) = data_included.Cluster_celllayer( data_included.Cluster_celllayer>2)+1; %>=3 -> n + 1 
%         data_included.Cluster_celllayer(data_included.Cluster_celldepth<split_4cab& data_included.Cluster_celllayer == 4) = 3; % lay 4cb
%     end
%     
%     % put WM before 6
%     temp = data_included.Cluster_celllayer;
%     wm_idx = max(temp);
%     temp(temp==wm_idx) = 0;
%     data_included.Cluster_celllayer = temp + 1;
%     
%     % flip array s.t. 6 is on bottom 
%     layers = unique(data_included.Cluster_celllayer);
%     layers = layers(~isnan(layers));
%     for i = 1:floor(length(layers)/2)
%         temp = data_included.Cluster_celllayer == i; 
%         data_included.Cluster_celllayer(data_included.Cluster_celllayer == length(layers)-i+1) = i;
%         data_included.Cluster_celllayer(temp) = length(layers)-i+1;
%     end
    