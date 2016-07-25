function [mu,covMat,params_updated] = updateCovMat(likelihood,givenData,regularization_param)

if nargin < 3
    regularization_param = 0;
end

params_updated = 1;
if sum(likelihood) == 0;
    mu = []; covMat = [];
    params_updated = 0;
    return;
end
    
  
dimension = size(givenData,2);

mu = sum(repmat(likelihood,1,dimension).*givenData)/sum(likelihood);

meanAdjustedData = givenData-repmat(mu,size(givenData,1),1);
meanAdjustedData = repmat(likelihood.^0.5,1,dimension).*meanAdjustedData;
covMat = (meanAdjustedData'*meanAdjustedData)/sum(likelihood) + eye(dimension)*regularization_param;

if ~isempty(find(isnan(covMat)))
    1
end

