% METHOD TO CLUSTER THE DATA GIVEN A GMC OBJECT
        % Inputs: obj = GMC object
        %         u   = N x d data to be clustered
        % Output: idx = N x 1 vector of cluster indices.
        function [idx,posteriors] = cluster(obj,u)
            
%             K = size(obj.mu,1);
%             d = size(obj.mu,2);
%             N = size(u,1);
            
            % Obtaining the inverse values with respect to the GMM marginal
            % distributions
            %inverseVals = computeInverseVals_vectorized(obj.mu,obj.Sigma,obj.PComponents,u,d,K,N,0);
            inverseVals = obj.computeInverseVals_vectorized(u,0);
            
            % Defining  the gmm object from which the gmc distribution is
            % derived
            gmmObject = gmdistribution(obj.mu,obj.Sigma,obj.PComponents);
            
            % Cluserting the inverse values using the gmm object.
            [idx,~,posteriors] = cluster(gmmObject,inverseVals);
        end
        