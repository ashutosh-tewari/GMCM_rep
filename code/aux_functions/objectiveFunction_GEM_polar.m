% This version of the objective function takes the inverse Cholesky factor as the input..
function [fVal,grad] = objectiveFunction_GEM_polar(x_in,sorting_order,u,d,M,N,post,g) %#codegen

%%%INPUT: x_in = the parameter vector
%%%Output: fVal = objectve function vlaue
% M = number of Modes
% d = dimension
% N = Number of data points
% u = NXd cdf value matrix (This remains fixed in the entire optimization process)


% GETTING THE MU, SIGMA AND ALPHA FROM THE PARAMTER VECTOR X_IN
% Getting the paramter values related to the means of modes from x_in
% vector

if nargin < 8
    g = 1000;
end

mean_params = [x_in(1:(M-1)*d)];
mean_param_mat = reshape(mean_params,d,M-1)';
mu(1,:) = zeros(1,d); % mu(2,:) = zeros(1,d); mu(2,1) = 1;
for j = 1:M-1
   radius = mean_param_mat(j,1);
   angles = mean_param_mat(j,2:end);
   mu(j+1,:) = polar2Cartesian(radius,angles)'; 
end
offSet1 = (M-1)*d;

% mean_params = [1 x_in(1:(M-1)*d-1)];
% mean_param_mat = reshape(mean_params,d,M-1)';
% mu = zeros(1,d);
% for j = 1:M-1
%    radius = mean_param_mat(j,1);
%    angles = mean_param_mat(j,2:end);
%    mu(j+1,:) = polar2Cartesian(radius,angles)'; 
% end
% % Removing these entries from the x_in vector
% offSet1 = (M-1)*d-1;
 
% mean_params = x_in(1:M*d);
% Reshaping it into a (KXd) matrix  
% mu = reshape(mean_params,M,d);
% Removing these entries from the x_in vector
% offSet1 = M*d;

% Getting the paramter values related to the choleskly factors of modes
num_W_params = (d+1)*d/2; % Number of params related to the cholesky factors 

W_mat = zeros(d,d,M);
Sigma = zeros(d,d,M);

temp_vec = cumsum(d:-1:1);

for k = 1:M
    % Getting choleskly factor parameters for the ith mode in the for loop.
    % Note: The first d elements in the vector corresponds to the diagonal
    % element and the remaining d*(d-1)/2 are the non-diagonal elements.
    current_mode_params = x_in(offSet1+(k-1)*num_W_params+1:offSet1+k*num_W_params);
    
    for j = 1:d
        W_mat(j:end,j,k) = current_mode_params(temp_vec(j)-d+j:temp_vec(j));
    end
    
    % Getting the covaraince matrix from a Cholskey factor
    Sigma(:,:,k) = inv(W_mat(:,:,k)'*W_mat(:,:,k)+ 1E-8*eye(d)); 
    
    clear current_mode_params diag_elements;
end

% % Removing the cholesky factor related entriers from the x_in vector
offSet2 = offSet1 + M*num_W_params; 

% Getting the Mixing weights.
%PComponents = x_in(offSet2+1:end);
PComponents = [x_in(offSet2+1:end) 1-sum(x_in(offSet2+1:end))];

% Obtaining the inverseValues
if isempty(sorting_order)
    z = computeInverseVals_vectorized(mu,Sigma,PComponents,u,d,M,N,0,g);
else
    z = computeInverseVals_vectorized(mu,Sigma,PComponents,u,d,M,N,1,g);
    for i = 1:d
       z(sorting_order(:,i),i) = z(:,i);
    end
end

% EVALUATING THE OBJECTIVE FUNCTION

p1 = repmat(log(PComponents),N,1);
Z_bar = repmat(z,[1 1 M]);
Z_bar = Z_bar-repmat(reshape(mu',[1 d M]),N,1);
p2=zeros(N,k);
p3 = zeros(N,k);
for k = 1:M
    p2(:,k) =  sum(log(diag(W_mat(:,:,k))));
    temp_mat = (W_mat(:,:,k)*Z_bar(:,:,k)')';
    p3(:,k) = sum(temp_mat.*temp_mat,2);
    clear temp_mat;
end
Q1 = sum(sum((p1+p2-p3/2).*post));

temp_mat = zeros(size(Z_bar));
for j = 1:d
    for k = 1:M
        temp_mat(:,j,k) = (PComponents(k)/sqrt(2*pi*Sigma(j,j,k)))*exp((-0.5/Sigma(j,j,k))*Z_bar(:,j,k).^2);
    end
end
Q2 = sum(sum(log(sum(temp_mat,3))));

% part3 = d*log(2*pi)/2;

fVal = -1*(Q1-Q2); %negative log likelihood

if nargout > 1  % provide gradient only when needed (does not compute gradient during function call for line search)
    % Gradient Computation
    interimParam.mu = mu;
    interimParam.Sigma = Sigma;
    interimParam.alpha = PComponents;
    interimParam.W = W_mat;
    interimParam.idx_tril = obtainIndexTriL(d);
    interimParam.idx_triu = setdiff(1:d^2,interimParam.idx_tril);
    interimParam.const_mat = sparse((eye(d^2) + obtainPermutationMatrix(d,d)));
    idx = sub2ind([d,d],1:d,1:d);
    A_idx = reshape(1:d*N,[d,N]);
    idx1 = kron(ones(d,1),A_idx);
    idx2 = kron(A_idx,ones(d,1));
    Z_bar_transposed = zeros(d,N,M);
    Wk_times_Z_bar_transposed_times = zeros(d,N,M);
    Zbar_kron_Wkzbar = zeros(d^2,N,M);
    for m=1:M
        interimParam.Precision(:,:,m) = interimParam.W(:,:,m)'*interimParam.W(:,:,m);
        interimParam.V(:,:,m) = inv(interimParam.W(:,:,m));
        mat1 = kron(eye(d),interimParam.W(:,:,m))*interimParam.const_mat;
        mat2 = -kron(interimParam.Sigma(:,:,m),interimParam.Sigma(:,:,m)');
        interimParam.dSig_rr_dW(:,:,m) =(mat1*mat2(:,idx));
        
        Z_temp = Z_bar(:,:,m)';
        Wk_Z_temp = interimParam.W(:,:,m)*Z_temp;
        Z_bar_transposed(:,:,m) = Z_temp;
        Wk_times_Z_bar_transposed_times(:,:,m) = Wk_Z_temp;        
        Zbar_kron_Wkzbar(:,:,m) = Z_temp(idx2).*Wk_Z_temp(idx1);
        
    end
   
    fracSample2Use = 1;
    selectedPoints = sort(randsample(N,N*fracSample2Use));
    % Initializing the gradient vectors    
    der_mu=zeros(d,M);der_W=zeros(d^2,M);
    for ii = 1:length(selectedPoints)
        idx = selectedPoints(ii);
        dz = getPartial_z(interimParam,z(idx,:)');
%         z_bar = reshape(Z_bar(ii,:,:),[d,M]);
        z_bar = reshape(Z_bar_transposed(:,idx,:),[d,M]);
        zbar_kron_wzbar = reshape(Zbar_kron_Wkzbar(:,idx,:),[d^2,M]);
        [~,dQ1_dmu,dQ1_dW] = getPartialDerQ1(interimParam,z_bar,dz,zbar_kron_wzbar,post(idx,:));
        [~,dQ2_dmu,dQ2_dW] = getPartialDerQ2_vectorized(interimParam,z_bar,dz);
%        der_alpha = der_alpha+dQ1_dalpha-dQ2_dalpha;
        der_mu = der_mu+dQ1_dmu-dQ2_dmu;
        der_W = der_W+dQ1_dW-dQ2_dW;
    end
    

    deriv_mu1 = reshape(reshape(der_mu,d,M)',M*d,1);
%    deriv_alpha1 = reshape(der_alpha,M,1);
    der_W = der_W(interimParam.idx_tril,:,:);
    deriv_W1 = reshape(der_W,M*d*(d+1)/2,1);

    % Getting the partial derivative w.r.t to alpha via finite difference
    delta = 1E-8;
    deriv_alpha1 = nan(M-1,1);
    for i = 1:M-1
        idx = M*d + M*d*(d+1)/2+i;
        x_plus = x_in;x_minus = x_in;
        x_plus(idx) = x_plus(idx) + delta;
        x_minus(idx) = x_minus(idx) - delta;
        val_plus = objectiveFunction_GEM(x_plus,sorting_order,u,d,M,N,post);
        val_minus = objectiveFunction_GEM(x_minus,sorting_order,u,d,M,N,post);
        deriv_alpha1(i) = (val_plus-val_minus)/(2*delta);    
    end
    deriv_alpha1 = -1*deriv_alpha1;
    
    % Concatenating the derivative to get the grad vector.
    grad = [deriv_mu1/fracSample2Use ; deriv_W1/fracSample2Use ;deriv_alpha1(1:M-1)]; % Division by fraction of sample used to scale up the gradient
    grad = (-1*grad);  % Negative sign due to negative log-likelihood.
    
end






