function [dQ2_dalpha,dQ2_dmu,dQ2_dW] = getPartialDerQ2_vectorized(param,zbar,dz)

% Input:  dz = struct with three fields (alpha= 1xdxM, mu=dxdxM, W=d^2xdxM
%         zbar = dx1 of inverse cdf values
% Output: dQ2_dalpha = 1xM vector
%         dQ2_dmu = dxM vector
%         dQ2_dW = d^2xM vector


[M,d] = size(param.mu);

% temp_mat dxM matrix
temp_mat = param.alpha(ones(1,d),:)/sqrt(2*pi);
for m = 1:M
   temp_mat(:,m) = temp_mat(:,m)./(diag(param.Sigma(:,:,m)).^0.5).*exp(-0.5*(zbar(:,m).^2)./diag(param.Sigma(:,:,m)));
end
F = sum(temp_mat,2);

temp_mat1 = zeros(d,M);
for m = 1:M
   temp_mat1(:,m) = temp_mat(:,m).*(zbar(:,m)./diag(param.Sigma(:,:,m)));
end
temp_vec1 = sum(temp_mat1,2);

% Getting the partial derivatives w.r.t to the parameters of each component
%% vectorized version


% M1 = (temp_mat./param.alpha(ones(1,d),:))./F(:,ones(1,M)) - reshape(dz.alpha,[d,M]).*V1(:,ones(1,M));
% dQ2_dalpha = sum(M1);
% clear M1;
dQ2_dalpha = zeros(1,M);

V1 = (temp_vec1./F);
M1 = temp_mat1./F(:,ones(1,M));
M2 = zeros(d,d,M);
for m = 1:M
    M2(:,:,m) = diag(M1(:,m));
end
M3 = -dz.mu.*V1(:,ones(1,d),ones(1,M)) + M2;
dQ2_dmu = reshape(sum(M3,2),[d,M]);
clear M1 M2 M3;

M1 = zeros(d,M);
for m = 1:M
    M1(:,m) = diag(param.Sigma(:,:,m));
end
M2 = (0.5*(temp_mat./(F(:,ones(1,M)).*M1)));
M3 = (0.5*(temp_mat1./(F(:,ones(1,M)).*M1))).*zbar;
M4 = reshape(M2,[1,d,M]);
M5 = reshape(M3,[1,d,M]);
V1p = V1';
M6 = -M4(ones(1,d^2),:,:).*param.dSig_rr_dW + M5(ones(1,d^2),:,:).*param.dSig_rr_dW - V1p(ones(1,d^2),:,ones(1,M)).*dz.W;
dQ2_dW = reshape(sum(M6,2),[d^2,M]);
clear M1 M2 M3 M4 M5 M6;

% %%
% for k = 1:M
%     % Getting the partial derivative w.r.t to alpha
%     for r = 1:d
%         dQ2_dalpha(k) = dQ2_dalpha(k) + (temp_mat(r,k)/param.alpha(k))/F(r) - (temp_vec1(r)/F(r))*dz.alpha(:,r,k);
%     end    
%     
%     % Getting the partial derivative w.r.t to mu
%     for r = 1:d
%         tv = zeros(d,1);tv(r) =  temp_mat1(r,k)/F(r);
%         dQ2_dmu(:,k) = dQ2_dmu(:,k) - (temp_vec1(r)/F(r))*dz.mu(:,r,k) + tv;
%     end  
%     
%     % Getting the partial derivative w.r.t to W
%     for r = 1:d
%         dQ2_dW(:,k) = dQ2_dW(:,k) - (temp_mat(r,k)/(2*F(r)*param.Sigma(r,r,k)))*param.dSig_rr_dW(:,r,k) + (temp_mat1(r,k)/(2*F(r)*param.Sigma(r,r,k)))*zbar(r,k)*param.dSig_rr_dW(:,r,k) - (temp_vec1(r)/F(r))*dz.W(:,r,k);
%     end 
% end


