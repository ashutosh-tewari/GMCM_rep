% METHOD TO COMPUTE THE PDF VALUES W.R.T THE GMC DISCTRIBUTION
% OBJECT
%Inputs:  obj = GMC object
%           u = N x d data to be clustered
%Output: pdf_vals = N x 1 vector of the pdf values.
function pdf_vals = pdf(obj,u)

    [K,d] = size(obj.mu);
    N = size(u,1);

    % Obtaining the inverse values with respect to the GMM marginal
    % distributions
%     inverseVals = computeInverseVals_vectorized(obj.mu,obj.Sigma,obj.PComponents,u,d,K,N,0);
    inverseVals = obj.computeInverseVals_vectorized(u,0);

    % Obtaining the log-likelihood of the numerator
    small_mat = 1e-323*ones(N,1);
    first_part = zeros(size(inverseVals,1),K);
    inverseVals_hat = zeros(N,d,K);
    for k = 1:K
        inverseVals_hat(:,:,k) = inverseVals - repmat(obj.mu(k,:),N,1); % Getting the mean adjusted inverse vals
        temp_mat = (obj.Sigma(:,:,k)\inverseVals_hat(:,:,k)')';
        first_part(:,k) = obj.PComponents(k)*(2*pi)^(-d/2)*(det(obj.Sigma(:,:,k)))^-0.5*exp(-0.5*sum(inverseVals_hat(:,:,k).*temp_mat,2));
        clear temp_mat;
    end
    first_part_ll = log(sum(first_part,2) + small_mat);  % A small positive number is added to avoid log(0)
    clear inverseVals;

    % Getting the log-likelihood of the denominator
    second_part = zeros(N,K);
    for j = 1:d
        temp_vector = zeros(N,K);
        for k = 1:K
            temp_vector(:,k) =  obj.PComponents(k)*(1/sqrt(2*pi*obj.Sigma(j,j,k)))...
                *exp(-0.5*(1/obj.Sigma(j,j,k))*(inverseVals_hat(:,j,k).^2));
        end
        second_part(:,j) = log(sum(temp_vector,2)+ small_mat);
    end
    second_part_ll = sum(second_part,2);
    clear inverseVals_hat;

    log_likelihood = first_part_ll - second_part_ll;
    pdf_vals = exp(log_likelihood);