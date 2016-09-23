function inverseVals = computeInverseValsGMM(gmmObj,u)

bin_num = 1000;


min_val = min(gmmObj.mu)-8*(max(gmmObj.Sigma)^0.5);
max_val = max(gmmObj.mu)+8*(max(gmmObj.Sigma)^0.5);
range_vec = (min_val:(max_val-min_val)/(bin_num-1):max_val)';

cdf_vals = cdf(gmmObj,range_vec);
inverseVals = interp1q(cdf_vals,range_vec,u);
