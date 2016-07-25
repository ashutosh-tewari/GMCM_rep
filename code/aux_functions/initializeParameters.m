function [x_init,bounds,A,b] = initializeParameters(obj,u,d,K,initOpt)


if isstruct(initOpt)
    mu = initOpt.mu;
    Sigma = initOpt.Sigma;
    PComponents = initOpt.PComponents;
elseif strcmp(initOpt,'randSample')
    mu = randn(K,d);
    Sigma = repmat(eye(d),[1 1 K]);
    PComponents = rand(1,K);PComponents = PComponents/sum(PComponents);
elseif strcmp(initOpt,'psuedoEM')
    numRep=1;
    numIterations = 500;
    [mu,Sigma,PComponents] = obj.psuedo_EM(u,K,numIterations,numRep,1);
else
   error('Invalid initialization option.'); 
end

mu_matrix = mu;
x_init = mu_matrix(1:K*d);

for k = 1:K
    for j = 1:d 
        inv_chol_mat = inv(chol(Sigma(:,:,k))');
        x_init = [x_init inv_chol_mat(j:d,j)'];
    end
end

x_init = [x_init PComponents(1:K-1)];


% Assigning Bounds to the means
lb = [];ub=[];
for i = 1:K
    lb = [lb;-inf(d,1)];
    ub = [ub;+inf(d,1)];
end

% Assigning Bounds to the Covariance matrix elements
diag_idx = 1:d+1:d^2;
lb_mat = -inf(d);
lb_mat(diag_idx) = 0.1;
lb_mat = tril(lb_mat);

% Defining the lower bounds
for k = 1:K
    for j = 1:d 
       lb= [lb;lb_mat(j:d,j)];
    end
end

% Defining the upper bounds
ub_mat = tril(inf(d));
for k = 1:K
    for j = 1:d 
       ub= [ub;ub_mat(j:d,j)];
    end
end

% Assigning Bounds to the mixing weights
lb = [lb;0.001*ones(K-1,1)];
ub = [ub;0.999*ones(K-1,1)];

bounds = [lb ub];


% Setting the constraint that that sum of PComponents should be less than unity.
A = [zeros(1,numel(x_init)-K+1) ones(1,K-1)];
b = 0.99;
 
