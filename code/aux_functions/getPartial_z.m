function [dz] = getPartial_z(param,z)

% Input:  dz = struct with three fields (alpha= 1xdxM, mu=dxdxM, W=d^2xdxM
%         z = dx1 of inverse cdf values
% Output: dz = struct with three fields (alpha= 1xdxM, mu=dxdxM, W=d^2xdxM
[M,d] = size(param.mu);

z_standardized = zeros(d,M);
for m = 1:M
    z_standardized(:,m) = (z - param.mu(m,:)')./(diag(param.Sigma(:,:,m)).^0.5);
end

temp_mat = param.alpha(ones(1,d),:)/sqrt(2*pi);
for m = 1:M
   temp_mat(:,m) = temp_mat(:,m)./(diag(param.Sigma(:,:,m)).^0.5);
end

% Denominator of the triple product (dx1 vector)
dTau0_dz = sum(temp_mat.*exp(-0.5*z_standardized.^2),2);
dTau0_dz = dTau0_dz';


% % Derivative w.r.t alpha
% dz.alpha = reshape(0.5*(1+ erf(z_standardized/sqrt(2))),[1,d,M]);

% Derivative w.r.t  mu
dz.mu = zeros(d,d,M);
for m = 1:M
    dz.mu(:,:,m) = diag(-temp_mat(:,m).*exp(-0.5*z_standardized(:,m).^2));
end

% Derivative w.r.t W
dz.W = zeros(d^2,d,M);
for m = 1:M
    temp_vec = -0.5*temp_mat(:,m).*exp(-0.5*z_standardized(:,m).^2).*(z-param.mu(m,:)');
    temp_vec = (temp_vec./diag(param.Sigma(:,:,m)))';
    dz.W(:,:,m) = param.dSig_rr_dW(:,:,m).*temp_vec(ones(1,d^2),:);
end

% Dividing with the denominator (triple product rule)
% dz.alpha = -dz.alpha./dTau0_dz(:,:,ones(1,M));
dz.mu = -dz.mu./dTau0_dz(ones(1,d),:,ones(1,M));
dz.W = -dz.W./dTau0_dz(ones(1,d^2),:,ones(1,M));

%Forcing the gradient corresponding to the upper triangular part to zero
dz.W(param.idx_triu,:,:)=0;







