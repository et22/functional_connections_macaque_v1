function [h,p, chi2stat,df] = prop_test(X , N)
    df=1; % 2 samples
    
    % Observed data
    n1 = X(1);
    n2 = X(2);
    N1 = N(1);
    N2 = N(2);
    
    % Pooled estimate of proportion
    p0 = (n1+n2) / (N1+N2);
    
    % Expected counts under H0 (null hypothesis)
    n10 = N1 * p0;
    n20 = N2 * p0;
    observed = [n1 N1-n1 n2 N2-n2];
    expected = [n10 N1-n10 n20 N2-n20];
    
    % Standard Chi-square test
    chi2stat = sum((observed-expected).^2 ./ expected);
    p = 1 - chi2cdf(chi2stat,1);

    h=0;
    if p<0.05
        h=1;
    end
end