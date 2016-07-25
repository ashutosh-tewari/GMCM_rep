% METHOD TO OBTAIN THE MARGINL GMC DISTRIBUTION SPECIFIED BY THE INDICES
%Inputs:  obj_out = GMC object
%           idx = a vector of indices on which the marginal is desired
%Output: obj_marg = marginal gmc density.
function obj_marg = obtainMarginal(obj,idx)

    mu = obj.mu;
    Sigma = obj.Sigma;
    PComponents = obj.PComponents;
    
    mu = mu(:,idx);
    Sigma = Sigma(idx,idx,:);
    
    obj_marg = gmcdistribution(mu,Sigma,PComponents);

end