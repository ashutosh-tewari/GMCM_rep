function [dQ1_dalpha,dQ1_dmu,dQ1_dW] = getPartialDerQ1(param,zbar,dz,kronMat,post)

% Input:  dz = struct with three fields (alpha= 1xdxM, mu=dxdxM, W=d^2xdxM
%         zbar = dx1 of inverse cdf values
% Output: dQ1_dalpha = 1xM vector
%         dQ1_dmu = dxM vector
%         dQ1_dW = d^2xM vector


[M,d] = size(param.mu);

% temp_mat will be used as an intermidate computation in the calculation of
% all the three partial derivatives.
temp_mat = zeros(d,1);
for l = 1:M
    temp_mat = temp_mat+ param.Precision(:,:,l)*zbar(:,l)*post(l);
end


% Getting the partial derivative w.r.t alpha, mu and W
dQ1_dalpha = zeros(1,M);
dQ1_dmu = zeros(d,M);
dQ1_dW = zeros(d^2,M); 
for k = 1:M
%     dQ1_dalpha(k) = post(k)/param.alpha(k) - dz.alpha(:,:,k)*temp_mat;
    dQ1_dmu(:,k) = -1* dz.mu(:,:,k)*temp_mat + param.Precision(:,:,k)*(zbar(:,k)*post(k));
    invW_tran = param.V(:,:,k)';
    Wk_zbar = param.W(:,:,k)*zbar(:,k);
%     V1 = kron(zbar(:,k),Wk_zbar);
    V1 = kronMat(:,k);
    dQ1_dW(:,k) = invW_tran(:)*post(k) - V1*post(k) - dz.W(:,:,k)*temp_mat;
end


%Forcing the gradient corresponding to the upper triangular part to zero
dQ1_dW(param.idx_triu,:)=0;








