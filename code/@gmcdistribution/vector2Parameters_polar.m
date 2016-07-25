function obj = vector2Parameters_polar(obj,params)

d = obj.NDimensions;
K = obj.NComponents;

% Defining the GMM object using the parameter vector
mean_params = [params(1:(K-1)*d)];
mean_param_mat = reshape(mean_params,d,K-1)';
mu(1,:) = zeros(1,d); %mu(2,:) = zeros(1,d); mu(2,1) = radius;
for j = 1:K-1
   radius = mean_param_mat(j,1);
   angles = mean_param_mat(j,2:end);
   mu(j+1,:) = polar2Cartesian(radius,angles)'; 
end
offSet1 = (K-1)*d;
% Removing these entries from the x_in vector

% % Defining the GMM object using the parameter vector
% mean_params = [1 params(1:(K-1)*d-1)];
% mean_param_mat = reshape(mean_params,d,K-1)';
% mu = zeros(1,d);
% for j = 1:K-1
%    radius = mean_param_mat(j,1);
%    angles = mean_param_mat(j,2:end);
%    mu(j+1,:) = polar2Cartesian(radius,angles)'; 
% end
% % Removing these entries from the x_in vector
% offSet1 = (K-1)*d-1;

% mu_params = params(1:K*d);
% mu = reshape(mu_params,K,d);
% offSet1 = K*d;

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
