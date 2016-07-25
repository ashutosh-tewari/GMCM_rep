% This function implements the newly derived Generalized EM algorithm to estimate the parameters of a
% Gaussian Mixture Copula model
 
function [x_final,maxLLVal,output_ll] = GEM(obj,u,K,x_init,numIterations,A,b,bounds,plo)

if nargin < 9
    plo = 0;
end

if plo
    figure;
end

[N,d] = size(u);

% sorting the CDF values for faster interpolation
for i = 1:d
    [u_sorted(:,i),sort_idx(:,i)] = sort(u(:,i));
end

options = optimset('Display','off','MaxFunEvals',500000,'DiffMinChange',1e-10,'TolFun',1e-10,'TolCon',1e-10,'TolX',1e-10,'MaxIter',1,'Algorithm','interior-point','GradObj','off');

ll_temp_1 = nan(1,numIterations);
Q_t1_t = nan(1,numIterations);
x_final = nan(numIterations,length(x_init));

for count1 = 1:numIterations
    obj = obj.vector2Parameters(x_init);
    ll_temp_1(count1) = sum(log(obj.pdf(u)));
    [~,post] = obj.cluster(u);

    % Updating the parameters 
%    Q_t_t(count1) = objectiveFunction_GEM(x_init,sort_idx,u_sorted,d,K,N,post); 
    [x_final(count1,:),Q_t1_t(count1)] = fmincon(@(x)objectiveFunction_GEM(x,sort_idx,u_sorted,d,K,N,post),x_init,A,b,[],[],bounds(:,1),bounds(:,2),[],options);
   
    subplot(4,1,1);plot(ll_temp_1,'ko-'); hold on;
    subplot(4,1,2);plot(Q_t1_t,'ro-');hold on;
    drawnow; 
   
    x_init = x_final(count1,:); 
end

x_final = x_init;
maxLLVal = sum(log(obj.pdf(u)));
