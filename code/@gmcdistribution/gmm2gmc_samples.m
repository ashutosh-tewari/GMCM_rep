%METHOD TO CONVERT THE SAMPLE OF UNDERLYING GMM TO THE GMC SAMPLES
        % Inputs: obj = GMC object
        %           N = number of samples required;
        % Output: samples = N x d matrix where d is the data dimension
        function gmc_samples = gmm2gmc_samples(obj,gmm_samples)
            
            K = size(obj.mu,1);
            d = size(obj.mu,2);
            
            % Defining  the gmm object from which the gmc distribution is
            % derived
            gmmObject = gmdistribution(obj.mu,obj.Sigma,obj.PComponents);
            
            % Obtaining the marginal distribution of the gmm.
            gmm_marginals = obtainMarginalsOfGMM(obj.mu,obj.Sigma,obj.PComponents,K,d);
            
            % samples from the gmc are nothing but the marginal cdf values
            % of the gmm samples.
            gmc_samples = nan(size(gmm_samples));
            for i=1:d
                gmc_samples(:,i) = cdf(gmm_marginals{i},gmm_samples(:,i));
            end
            
        end