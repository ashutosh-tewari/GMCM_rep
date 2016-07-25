% METHOD TO COMPUTE THE CDF VALUES W.R.T THE GMC DISCTRIBUTION
% OBJECT
%Inputs:  obj = GMC object
%           u = N x d data to be clustered
%Output: cdf_vals = N x 1 vector of the pdf values.
function cdf_vals = cdf(obj,u)

    K = size(obj.mu,1);
    d = size(obj.mu,2);
    N = size(u,1);

    % Obtaining the inverse values with respect to the GMM marginal
    % distributions
    inverseVals = computeInverseVals_vectorized(obj.mu,obj.Sigma,obj.PComponents,u,d,K,N,0);

    % Defining  the gmm object from which the gmc distribution is
    % derived
    gmmObject = gmdistribution(obj.mu,obj.Sigma,obj.PComponents);

    cdf_vals = cdf(gmmObject,inverseVals);

end