function postprocess()
clearvars; 
addpath(genpath(pwd));

flag = config();

for group_idx = 1:flag.group_cnt

% load and combine data from all sessions
num_neus = zeros(flag.num_sessions(group_idx),1);
num_pairs = zeros(flag.num_sessions(group_idx),1);
for i = 1:flag.num_sessions(group_idx)
    clab = "c" + int2str(i);
    session_number = i-1+flag.group_min_ses(group_idx);
    op.(clab) = load(flag.output_dir + flag.ccg_output_file + int2str(session_number)+".mat", 'ccg_output');
    % double check this with sessions 4 and 5
    op.(clab).ccg_output.data = split_layers_at_mean(op.(clab).ccg_output.data, flag, session_number); 
    op.(clab).ccg_output.data.Cluster_celldepth_scale = op.(clab).ccg_output.data.Cluster_celldepth;
    num_neus(i) = length(op.(clab).ccg_output.data.Cluster_celltype);
    num_pairs(i) = size(op.(clab).ccg_output.ccg_control,1);
end

% combining data across sessions
vert_attributes = ["ccg_control", "ccg_unnorm", "ccg_norm_jitter", "ccg_norm"];
vert_data_attributes = ["Cluster_MI_max", "Cluster_simpcomp"];
horiz_data_attributes = ["Cluster_celltype", "Cluster_celllayer", "Cluster_celldepth", "Cluster_celldepth_scale"];


for i = 2:flag.num_sessions(group_idx)
    clab = "c" + int2str(i);
    for ii = 1:length(vert_attributes)
        op.c1.ccg_output.(vert_attributes(ii)) = [op.c1.ccg_output.(vert_attributes(ii)); op.(clab).ccg_output.(vert_attributes(ii))];
    end
    for ii = 1:length(vert_data_attributes)
        op.c1.ccg_output.data.(vert_data_attributes(ii)) = [op.c1.ccg_output.data.(vert_data_attributes(ii)); op.(clab).ccg_output.data.(vert_data_attributes(ii))];
    end
    for ii = 1:length(horiz_data_attributes)
         op.c1.ccg_output.data.(horiz_data_attributes(ii)) = [op.c1.ccg_output.data.(horiz_data_attributes(ii)), op.(clab).ccg_output.data.(horiz_data_attributes(ii))];
    end
end

% add 1 to cell type, so we start counting from 0 instead of 1
op.c1.ccg_output.data.Cluster_celltype  = op.c1.ccg_output.data.Cluster_celltype+1;

op.c1.ccg_output.session = [];
op.c1.ccg_output.data.Cluster_session = [];

conditionID = {};
spiketrain = {};

for i = 1:flag.num_sessions(group_idx)
    clab = "c" + int2str(i);
    op.c1.ccg_output.session = [op.c1.ccg_output.session, i*ones(1,num_pairs(i))];
    op.c1.ccg_output.data.Cluster_session = [op.c1.ccg_output.data.Cluster_session , i*ones(1,num_neus(i))];
    
    conditionID = [conditionID; op.(clab).ccg_output.data.conditionID];
    spiketrain = [spiketrain; op.(clab).ccg_output.data.spiketrain];
end

op.c1.ccg_output.data.conditionID  = conditionID;
op.c1.ccg_output.data.spiketrain = spiketrain;

for i = 2:flag.num_sessions(group_idx)
    clab = "c" + int2str(i);
    op.c1.ccg_output.neuron_id_pairs = [op.c1.ccg_output.neuron_id_pairs; op.(clab).ccg_output.neuron_id_pairs + sum(num_neus(1:i-1))];
end

ccg_output = op.c1.ccg_output;
ccg_output.data.layer_centers = [1,2,3,4,5,6,7];

% get ccg attributes
ccg_data = ccg_to_attributes(ccg_output, flag);
ccg_data.compute_flag = op.c1.ccg_output.flag;
ccg_data.postproc_flag = flag;
save(flag.postproc_output(group_idx), 'ccg_data');
end