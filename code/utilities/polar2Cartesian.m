function x = polar2Cartesian(r,angles)
%UNTITLED3 Summary of this function goes here
%  function that converts the polar coordinates into cartesian coordinates

    % adding small perturbation to angles
    angles = angles+eps;
    
    dim = length(angles)+1;
    x = nan(dim,1);
    
         
    temp_mat1 = zeros(dim);
    temp_mat1 = temp_mat1 + diag([cos(angles) 0]);
    temp_mat1(~temp_mat1) = 1;
    
    temp_mat2 = nan(dim);
    temp_mat2(2:end,1:end-1) = repmat(sin(angles),dim-1,1);    
    temp_mat2 = tril(temp_mat2,-1);
    temp_mat2(~temp_mat2) = 1;
    
    temp_mat = temp_mat1.*temp_mat2;
    
    x = r*prod(temp_mat,2);


end

