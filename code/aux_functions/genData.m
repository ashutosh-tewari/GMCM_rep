function [data_out,mu1,Sigma,PComponents] = genData(N,d,K,plo)

if nargin < 4
    plo = 0;
end

for k = 1:K
   mu(k,:) = unidrnd(12,1,d);
   temp_mat = randn(d);
   Sigma(:,:,k) = temp_mat'*temp_mat;
   PComponents(k) = rand(1);
end
mu1 = mu;

PComponents =  PComponents/sum(PComponents);
gmmObj = gmdistribution(mu,Sigma,PComponents);

for j = 1:d
    marg_dist1{j} = gmdistribution(mu(:,j),Sigma(j,j,:),PComponents);
end

data = random(gmmObj,N);
clear mu sig PC;
for j = 1:d
   for k = 1:unidrnd(2)
       mu(k) = 2*unidrnd(20);
       sig(k) = unifrnd(1,4);
       PC(k) = rand(1);
   end
   PC = PC/sum(PC);
   marg_dist2{j} = gmdistribution(mu',reshape(sig,1,1,numel(PC)),PC); 
end

for j = 1:d
    u(:,j) = cdf(marg_dist1{j},data(:,j));
    data_out(:,j) = computeInverseValsUnivGMM(marg_dist2{j},u(:,j));
end

if plo
    figure;
    subplot(1,2,1);plot(u(:,1),u(:,2),'ko');
    subplot(1,2,2);plot(data_out(:,1),data_out(:,2),'ko');
end