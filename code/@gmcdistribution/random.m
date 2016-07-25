%METHOD TO SAMPLE FROM THE GMC DISTRIBUTION 
        % Inputs: obj = GMC object
        %           N = number of samples required;
        % Output: samples = N x d matrix where d is the data dimension
        function gmc_samples = random(obj,N)
            
            K = size(obj.mu,1);
            d = size(obj.mu,2);
            
            % Defining  the gmm object from which the gmc distribution is
            % derived
            gmmObject = gmdistribution(obj.mu,obj.Sigma,obj.PComponents);
            % Sampling from the gmm object
            gmm_samples = random(gmmObject,N);
            
            % Obtaining the marginal distribution of the gmm.
            gmm_marginals = obtainMarginalsOfGMM(obj.mu,obj.Sigma,obj.PComponents,K,d);
            
            % samples from the gmc are nothing but the marginal cdf values
            % of the gmm samples.
            gmc_samples = nan(size(gmm_samples));
            for i=1:d
                gmc_samples(:,i) = cdf(gmm_marginals{i},gmm_samples(:,i));
            end
            
        end
        
        
        % This function computes the univariate GMM object from the multivariate GMM params.
        function marginalsGMM = obtainMarginalsOfGMM(mu,Sigma,PComponents,K,d)

        for i = 1:d
            for j = 1:K
                mu_marginal(j,1)  = mu(j,i);
                sigma_marginal(:,:,j) = Sigma(i,i,j);
                alpha_marginal(j) = PComponents(j)+10^-300;
            end    
            marginalsGMM{i} = gmdistribution(mu_marginal,sigma_marginal,alpha_marginal);
            clear mu_marginal sigma_marginal alpha_marginal;
        end
        end
        