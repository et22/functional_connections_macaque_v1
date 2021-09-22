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