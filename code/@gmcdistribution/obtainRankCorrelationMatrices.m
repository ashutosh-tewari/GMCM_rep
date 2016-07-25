        % METHOD TO compute the Spearman Correlation Coefficients
        % Inputs: obj = GMC object
        %         u   = N x d data to be clustered
        %         rank_matrix = N x d matrix consisting of columnwise ranks of entries
        % Output: spearman_coeff = d x d x K (matrices of spearman correlation coefficient).
        function [rankCorrelationMatrix,new_obj] = obtainRankCorrelationMatrices(obj,u)
            
            K = size(obj.mu,1);
            d = size(obj.mu,2);
            N = size(u,1);
            
            rankCorrelationMatrix = nan(d,d,K); % defining the 3-D rank correlation matrix.
            
            % Obtaining the inverse values with respect to the GMM marginal
            % distributions
            inverseVals = computeInverseVals_vectorized(obj.mu,obj.sigma,obj.alpha,u,d,K,N,0);            
            
           for k = 1:K
               likelihood(:,k) = mvnpdf(inverseVals,obj.mu(k,:),obj.sigma(:,:,k));
           end
           likelihood = likelihood.*repmat(obj.alpha,N,1);
           likelihood =  (likelihood./repmat(sum(likelihood,2),1,K));

%            end   
           
            gmmObject = gmdistribution(obj.mu,obj.sigma,obj.alpha);            
            % Cluserting the inverse values using the gmm object.
            idx = cluster(gmmObject,inverseVals);
           for k = 1:K
               u_cluster = u(idx==k,:);
               for j = 1:d
                    [~,idx1] = sort(u_cluster(:,j));
                    [~,idx2] = sort(idx1);
                    rank_matrix_cluster(:,j) = idx2;
                    clear idx1 idx2;
               end
               meanAdjustedRankCluster = rank_matrix_cluster - repmat(mean(rank_matrix_cluster),size(rank_matrix_cluster,1),1);
               temp1 = meanAdjustedRankCluster'*meanAdjustedRankCluster;
               diag_element = diag(temp1);
               temp2 = (diag_element*diag_element').^0.5;
               rankCorrelationMatrix(:,:,k) = temp1./temp2;
               clear rank_matrix_cluster  meanAdjustedRankCluster u_cluster temp1 temp2 diag_element;
           end
           
           rankCorrelationMatrix = 2*sin(rankCorrelationMatrix*(pi/6));  % Conversion from spearman rho to pearson correlation coeff. 
           
           for k = 1:K
               sigma_updated(:,:,k) = ((diag(obj.sigma(:,:,k))*diag(obj.sigma(:,:,k))').^0.5).*rankCorrelationMatrix(:,:,k);
           end
           
           new_obj = gmcdistribution(obj.mu,sigma_updated,obj.alpha);


        end
    