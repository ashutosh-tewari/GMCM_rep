% METHOD TO FIT A GMC DISTRIBUTION TO DATA
        % Inputs: u  :N x d matrix of CDF values obtained after fitting marginal densities.
        %       : K  : number of clusters
        %       : d  : number of dimensions
        %       : N  : number of data samples.
        % varargin:  'Start' :  'rand_init' OR 'EM_init'  (default'EM_init');
        %            'replicates' : Integer  (default 1)
        %            'iteration'  : Integer  (default 100)
        %            'algorithm'  : 'active-set' OR 'interior-point' (default 'active-set')
        % Example:
        % obj.fit(u,K,d,N,'Start','rand_init','replicates',20,'algorithm','interior-point')
        
        function [obj] = fit(obj,u,K,d,N,varargin)
           
            % parsing the input argument
            p = inputParser;
            p.addParamValue('Start','randSample');
            p.addParamValue('replicates',1);
            p.addParamValue('iteration_GD',0);
            p.addParamValue('iteration_EM',100);
            p.addParamValue('EM_rep',1);
            p.addParamValue('algorithm','active-set');
            p.addParamValue('grad_type','finite-difference');
            p.parse(varargin{:});
            
                        
            for i = 1:p.Results.EM_rep
                % Checking the user specified initialization.
                [x_init,bounds,A,b] = initializeParameters(obj,u,d,K,p.Results.Start);
                % Implementing the Generalized EM algorithm
%                 [x_GEM_mat(:,i),maxLLVal(i,1)] = gmcm_GEM(u,K,d,N,x_init,p.Results.iteration_EM,A,b,bounds,1);
                [x_GEM_mat(:,i),maxLLVal(i,1)] = obj.GEM(u,K,x_init,p.Results.iteration_EM,A,b,bounds,1);
            end
            
            [~,max_idx] = max(maxLLVal);
            x_final = x_GEM_mat(:,max_idx);


            % converting the paramters in vectorized form into arrays
%             [mu,Sigma,PComponents] = vector2GMMParameters(x_final',d,K);
            obj = obj.vector2Parameters(x_final');

            
%             % updating the propoerties of the gmcdistribution object
%             obj.mu = mu;
%             obj.Sigma = Sigma;
%             obj.PComponents = PComponents;     
            
            obj.nLogL = -sum(log(obj.pdf(u)));
            obj.N = N;
            obj.numParams = length(x_final);
                        
        end