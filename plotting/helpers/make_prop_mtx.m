function [prop_mtx, prop_numer, prop_denom] = make_prop_mtx(cat_row, cat_col, numer_var, flip)
if ~exist('flip', 'var')
    flip = 0 ;
end
    cat_rows = unique(cat_row);
    cat_cols = unique(cat_col);
    cat_rows = cat_rows(~isnan(cat_rows));
    cat_cols = cat_cols(~isnan(cat_cols));
    for i = 1:length(cat_rows)
        for j = 1:length(cat_cols)
            prop_numer(i,j) = sum(numer_var(cat_row == cat_rows(i) & cat_col == cat_cols(j)));
            if flip && i~= j
                prop_denom(i,j) = numel(numer_var(cat_row == cat_rows(i) & cat_col == cat_cols(j))) + numel(numer_var(cat_row == cat_cols(i) & cat_col == cat_rows(j)));
            else
                prop_denom(i,j) = numel(numer_var(cat_row == cat_rows(i) & cat_col == cat_cols(j)));
            end
        end
    end
    prop_mtx = prop_numer./prop_denom;
end