% figuring out significance error
flag = config();
group_idx = 1;
load(flag.postproc_output(group_idx), 'ccg_data');

% subsetting ccg_data based on significance
ccg_data = ccg_data.ccg;
[ccg_data, sig_idx] = get_significant_ccgs(ccg_data, flag);

