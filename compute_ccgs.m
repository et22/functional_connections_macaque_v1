function compute_ccgs(ses_num)
addpath(genpath(pwd));

% get input configuration
flag = config();

% make output directory
if ~exist(flag.output_dir, 'dir')
    mkdir(flag.output_dir);
end

% load data
load(flag.data_dir + "ori1_data_session_" + int2str(ses_num) +".mat", 'data');
load(flag.data_dir + "ori1_MI_session " + int2str(ses_num) + ".mat", 'MI_max');

% subset data, only include visually responsive neurons
data_included.conditionID = data.conditionID;
data_included.spiketrain = data.spiketrain(:,:,data.Cluster_Included);
data_included.Cluster_celltype = data.Cluster_celltype(data.Cluster_Included);
data_included.Cluster_celldepth = data.Cluster_celldepth(data.Cluster_Included);
data_included.Cluster_celllayer = data.Cluster_celllayer(data.Cluster_Included);
data_included.Cluster_MI_max = MI_max(data.Cluster_Included);
data_included.Cluster_simpcomp = (data_included.Cluster_MI_max>1)+1; % 1 is complex, 2 is simple

if flag.debug_mode
    data_included = get_debug_subset(data_included);
end

ccg_output = analyze_ccgs(data_included, flag);
save(flag.output_dir + flag.ccg_output_file + int2str(ses_num)+".mat", 'ccg_output');
end

% helper function for debugging
function data_included = get_debug_subset(data_included)
    data_included.conditionID = data_included.conditionID(1:100);
    data_included.spiketrain = data_included.spiketrain(1:100,:,1:6);
    data_included.spiketrain = data_included.spiketrain(:,:,1:6);
    data_included.Cluster_celltype = data_included.Cluster_celltype(1:6);
    data_included.Cluster_celldepth = data_included.Cluster_celldepth(1:6);
    data_included.Cluster_celllayer = data_included.Cluster_celllayer(1:6);
    data_included.Cluster_MI_max = data_included.Cluster_MI_max(1:6);
    data_included.Cluster_simpcomp = data_included.Cluster_simpcomp(1:6); % 1 is complex, 2 is simple
end
