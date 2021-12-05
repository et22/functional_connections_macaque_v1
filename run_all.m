% compute ccgs and run all analyses
flag = config();

% computing
for i = 1:flag.num_sessions
    compute_ccgs(i)
end

% extracting attributes
postprocess()

%clustering
compute_clusters()

%plotting
plot_figure4()
plot_figure5()
plot_figures67()