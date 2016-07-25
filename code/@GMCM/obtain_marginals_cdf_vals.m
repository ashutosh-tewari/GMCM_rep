function u_mat = obtain_marginals_cdf_vals(obj,givenData) 

[N,d] = size(givenData);
u_mat = nan(N,d);

for i = 1:d
    u_mat(:,i) = obj.marginals{i}.cdf(givenData(:,i));
end
