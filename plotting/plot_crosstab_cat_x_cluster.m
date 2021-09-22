function plot_crosstab_cat_x_cluster(pre_cat, post_cat, cluster_idx)
% remove nans from categorical variables before running this fxn, or label
% them something else
cat_pairs = sortrows([pre_cat, post_cat]);
cat_pairs_idx = cantor_pairing_2to1(cat_pairs(:,1), cat_pairs(:,2));

[tbl,chi2,p,labels] = crosstab(cat_pairs_idx, cluster_idx);
disp("crosstab chi^2 = " + chi2 + ", p = " + p);