function obj = vector2Parameters(obj,params)

d = obj.NDimensions;
K = obj.NComponents;

% Defining the GMM object using the parameter vector
mu_params = params(1:K*d);
mu = reshape(mu_params,K,d);
offSet1 = K*d;

temp_vec = cumsum(d:-1:1);

num_W_params = (d+1)*d/2;
Sigma = zeros(d,d,K);
for k = 1:K
    current_mode_params = params(offSet1 + (k-1)*num_W_params+1 : offSet1 + k*num_W_params);
    W_mat = zeros(d,d);
    for j = 1:d
        W_mat(j:end,j) = current_mode_params(temp_vec(j)-d+j:temp_vec(j));
    end
    Sigma(:,:,k) = inv(W_mat'*W_mat);
    clear current_mode_params;
end

offSet2 = offSet1 + K*num_W_params;

PComponents = [params(offSet2+1:end) 1-sum(params(offSet2+1:end))];

% assigning the values
obj.mu = mu;
obj.Sigma = Sigma;
obj.PComponents = PComponents;
