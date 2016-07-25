% Getting the inverse values of univariate gmm CDFs empirically

function inverseVals = computeInverseVals_vectorized(obj,u,is_u_sorted,bin_num)

[N,d] = size(u);
K = obj.NComponents;

if nargin < 4
    bin_num = 1000;
end

mu_mat = zeros(1,d,K);
sigma_mat = zeros(1,d,K);
alpha_mat = zeros(1,d,K);
range_mat = zeros(bin_num,d);


for j = 1:d
    mu_mat(1,j,:) = reshape(obj.mu(:,j),1,1,K);
    sigma_mat(1,j,:) = (2*obj.Sigma(j,j,:)).^0.5;
    alpha_mat(1,j,:) = reshape(obj.PComponents,1,1,K)/2;
    min_val = floor(min(obj.mu(:,j))-8*(max(obj.Sigma(j,j,:))^0.5));
    max_val = ceil(max(obj.mu(:,j))+8*(max(obj.Sigma(j,j,:))^0.5));
    range_mat(:,j) = (min_val:(max_val-min_val)/(bin_num-1):max_val)';
end
mu_mat = reshape(obj.mu',1,d,K);

erf_val_x = erf((range_mat(:,:,ones(K,1)) - mu_mat(ones(bin_num,1),:,:))./sigma_mat(ones(bin_num,1),:,:));

uVals = sum(alpha_mat(ones(bin_num,1),:,:).*(ones(bin_num,d,K)+erf_val_x),3);

inverseVals = zeros(N,d);
for j = 1:d 
    if is_u_sorted
        inverseVals(:,j) = interp1q_custom(uVals(:,j),range_mat(:,j),u(:,j));
    else
        inverseVals(:,j) = interp1q(uVals(:,j),range_mat(:,j),u(:,j));
    end
end
