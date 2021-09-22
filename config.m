function flag = config()
%debugging?
flag.debug_mode = 0;

%number of sessions
flag.num_sessions = 5;

%separate or group sessions
flag.group_sessions = false;
flag.group_by_animals = false;

if ~flag.group_sessions
    flag.group_cnt = 1;
    flag.group_min_ses = [1];
    flag.group_max_ses = [5];
    flag.group_labels = [""];
    flag.num_sessions = [5];
elseif flag.group_by_animals
    flag.group_cnt = 2;
    flag.group_min_ses = [1 4];
    flag.group_max_ses = [3 5];
    flag.group_labels = ["monkey 1", "monkey 2"];
    flag.num_sessions = [3,2];
else
    flag.group_cnt = 5;
    flag.group_min_ses = [1 2 3 4 5];
    flag.group_max_ses = [1 2 3 4 5];
    flag.num_sessions = [1,1,1,1,1];
    flag.group_labels = ["ses 1", "ses 2", "ses 3", "ses 4", "ses 5"];
end

%ccg calculation variables
flag.start_time = 500; % start 400 ms after stimulus onset
flag.end_time = 1100; % end 1000 ms after stimulus onset
flag.max_lag = 100;
flag.min_lag = -100;
flag.jit_window = 25;

%postprocess variables
flag.max_pt_lag = 10;
flag.min_pt_lag = -10;
flag.split_layers = 1;
flag.large_output = 1;

%sig criteria
flag.sig_num_stds = 7; % must be 7 stds above noise mean
flag.sig_max_lag = 10; % must have peak w/in 10 ms of 0
flag.sig_min_std = 0; % must have a non-zero noise std

%cluster variables
flag.rng_seed = 1;
flag.tsne_dims = 3;
flag.max_k = 10;
flag.kmeans_rep = 50;
flag.kmeans_maxiter = 100;

%where to save
flag.filename = "";
flag.date = '25-Aug-2021';%date;
if ~flag.debug_mode
    flag.output_dir = "output/" + flag.date + "/";
else
    flag.output_dir = "output/debug/";
end
flag.data_dir = "data/";
flag.ccg_output_file = "ccg_output_";

for i = 1:flag.group_cnt
    if flag.large_output
        flag.postproc_output(i) = flag.output_dir + "ccg_attributes" + "_large_" + flag.group_labels(i) + ".mat";
    else
        flag.postproc_output(i) = flag.output_dir + "ccg_attributes" + "_small_" + flag.group_labels(i) + ".mat";
    end
    flag.cluster_output(i) = flag.output_dir + "cluster_output" + flag.group_labels(i) + ".mat";
end
flag.figure_dir = "figures/matlab/";
