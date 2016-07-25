% Function to implement psuedo EM algorithm
function [mu_final,sigma_final,alpha_final] = psuedo_EM(obj,u,K,numIterations,numRep,plo)


if nargin < 6
    plo = 0;
end

if plo
    figure;
end

[N,d] = size(u);

for c = 1:numRep
    [x_init,~,~,~] = initializeParameters(obj,u,d,K,'randSample');

%     [obj.mu,obj.Sigma,obj.PComponents] = vector2GMMParameters(x_init,d,K);
    obj = vector2Parameters(obj,x_init);

    % Getting the inverse values of univariate gmm CDFs empirically for the
    % first iteration.
%     inverseVals = obj.computeInverseVals_vectorized(u,0);

    maxLLVal{c}  = -inf;
    for count = 1:numIterations   

%         marginalGMM = obtainMarginalsOfGMM(mu,Sigma,PComponents,K,d);
        %%%STEP 1
        try
            ll_temp = sum(log(obj.pdf(u)));
        catch
            ll_temp = -inf;
        end

        %ll_temp1 = obtainCopulaLikelihood_V1(mu,Sigma,PComponents,inverseVals,d,K,N);

        %%%STEP 2
        if ll_temp > maxLLVal{c}
            clear inverseVals;
            % updating the inverse values
            inverseVals = obj.computeInverseVals_vectorized(u,0);
            % updated theta final
            mu_rep{c} = obj.mu;
            sigma_rep{c} = obj.Sigma;
            alpha_rep{c} = obj.PComponents;
            %updated maxLLVal
            maxLLVal{c} = ll_temp;

        end

        if plo
            plot(count,maxLLVal{c},'r.-');hold on;
            drawnow;
        end 

        for j = 1:K        
             inverse_vals_hat = inverseVals - repmat(obj.mu(j,:),N,1); % Getting the mean adjusted inverse vals
             temp_mat = (obj.Sigma(:,:,j)\inverse_vals_hat')';%temp_mat = y_hat(:,:,k)*(inv(V_mat(:,:,k)))';
             temp_vec = -(d/2)*log(2*pi) -log(det(obj.Sigma(:,:,j))) - 0.5*sum(inverse_vals_hat.*temp_mat,2);
             likelihood(:,j) = exp(temp_vec);
             clear temp_mat temp_vec inverse_vals_hat
        end


        %%%STEP 3
        % E Step
        likelihood = likelihood.*repmat(obj.PComponents,size(likelihood,1),1) + 10^-300; % adding a small number just to avoid a all zero likelihood condition
        likelihood =  (likelihood./repmat(sum(likelihood,2),1,K));


        % M Step
        for j = 1:K
            [muUpdate,sigmaUpdate,params_updated] = updateCovMat(likelihood(:,j),inverseVals,0.001);
            if params_updated
                obj.mu(j,:) = muUpdate;
                obj.Sigma(:,:,j) = sigmaUpdate;
            end
            clear muUpdate sigmaUpdate;
        end

        PComponents = sum(likelihood)/size(u,1); 
        obj.PComponents = PComponents/sum(PComponents);

        clear likelihood;

        output_ll(count,1) = maxLLVal{c};

    end

end

[~,max_id] = max(cell2mat(maxLLVal));
mu_final = mu_rep{max_id};
sigma_final = sigma_rep{max_id};
alpha_final = alpha_rep{max_id};

