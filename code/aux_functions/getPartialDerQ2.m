function [dQ2_dalpha,dQ2_dmu,dQ2_dW] = getPartialDerQ2(param,zbar,dz)

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
dQ2_dalpha = zeros(1,M);
dQ2_dmu = zeros(d,M);
dQ2_dW = zeros(d^2,M);
for k = 1:M
%     % Getting the partial derivative w.r.t to alpha
%     for r = 1:d
%         dQ2_dalpha(k) = dQ2_dalpha(k) + (temp_mat(r,k)/param.alpha(k))/F(r) - (temp_vec1(r)/F(r))*dz.alpha(:,r,k);
%     end    
    
    % Getting the partial derivative w.r.t to mu
    for r = 1:d
        tv = zeros(d,1);tv(r) =  temp_mat1(r,k)/F(r);
        dQ2_dmu(:,k) = dQ2_dmu(:,k) - (temp_vec1(r)/F(r))*dz.mu(:,r,k) + tv;
    end  
    
    % Getting the partial derivative w.r.t to W
    for r = 1:d
        dQ2_dW(:,k) = dQ2_dW(:,k) - (temp_mat(r,k)/(2*F(r)*param.Sigma(r,r,k)))*param.dSig_rr_dW(:,r,k) + (temp_mat1(r,k)/(2*F(r)*param.Sigma(r,r,k)))*zbar(r,k)*param.dSig_rr_dW(:,r,k) - (temp_vec1(r)/F(r))*dz.W(:,r,k);
    end 
end


% %% Getting the partial derivative w.r.t alpha
% first_part = temp_mat./param.alpha(ones(1,d),:);
% first_part = first_part./F(:,ones(1,M));
% first_part = sum(first_part);
% 
% temp_mat1 = zeros(d,M);
% for m = 1:M
%    temp_mat1(:,m) = temp_mat(:,m).*(zbar(:,m)./diag(param.Sigma(:,:,m)));
% end
% 
% for m = 1:M    
%     second_part(m) = sum(sum(temp_mat1,2).*(dz.alpha(:,:,m)'));
% end
% 
% dQ2_dalpha = first_part - second_part;
% clear first_part second_part;
% 
% %% Getting partial derivative w.r.t mu
% first_part = zeros(d,M);
% for m = 1:M
%     dzbar_dmu = dz.mu(:,:,m) -eye(d);
%     first_part(:,m) = temp_mat(:,m)./diag(param.Sigma(:,:,m)).*diag(dzbar_dmu);
% end
% 
% dQ2_dmu = zeros(d,M);
% for m = 1:M
%     dQ2_dmu(:,m) = -(zbar(:,m)./F).*sum(first_part,2);
% end
% 
% %% Getting the derivative w.r.t W
% mat = zeros(d,M);
% for m = 1:M
%     mat = temp_mat(:,m)./diag(param.Sigma(:,:,m));
% end
% mat = sum(mat,2);
% 
% third_part = zeros(d^2,M);
% for m = 1:M
%     vec= zeros(d^2,M);
%     for j = 1:d
%         vec(:,j) = (zbar(j)/F(j))*dz.W(:,j,m)*mat(j);
%     end
%     third_part(:,m) = sum(vec,2);
% end
% 
% second_part = zeros(d^2,M);
% for m = 1:M
%    mat = zeros(d^2,M);
%    for j = 1:d
%       mat(:,j) = (temp_mat(j,m)*zbar(j)^2)/(2*F(j)*param.Sigma(j,j,m)^2)*param.dSig_rr_dW(:,j,m); 
%    end
%    second_part(:,m) = sum(mat,2);
% end
% 
% first_part = zeros(d^2,M);
% for m = 1:M
%    mat = zeros(d^2,M);
%    for j = 1:d
%       mat(:,j) = temp_mat(j,m)/(2*F(j)*param.Sigma(j,j,m))*param.dSig_rr_dW(:,j,m); 
%    end
%    first_part(:,m) = sum(mat,2);
% end
% 
% dQ2_dW = -first_part + second_part - third_part;
% 
