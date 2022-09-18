addpath(genpath(pwd));

% compute ccgs and run all analyses
flag = config();

% save config file to load in python
save(append(flag.output_dir, 'config.mat'), '-struct', 'flag');

% computing
for i = 1:flag.num_sessions
    compute_ccgs(i)
end

% extracting attributes
postprocess()

%clustering
compute_clusters()

%plotting
plot_all()

%all other plots can be generated using the plot_figues123.ipynb folder