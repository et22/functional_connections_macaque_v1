function p = z_prop_test(null_prop, actual_prop, n)
    z_stat = (actual_prop-null_prop)/sqrt((null_prop*(1-null_prop))/(n));
    p = 1 - normcdf(abs(z_stat), 0, 1) + normcdf(-abs(z_stat),0,1);
end