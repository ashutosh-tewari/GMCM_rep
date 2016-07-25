% Author: Ashutosh Tewari (tewaria@utrc.utc.com)
% Affiliation: Decision Support & Machine Intelligence
%              United Technologies Research Center
%              East Hartford, CT: 06118
% Code based on the paper 'Parametric Characterization of MultiModal Densities with non-Gaussian
% Nodes', A. Tewari, M.J. Giering, A. Raghunathan, OEDM workshop, ICDM 2011 (paper included in the package) 

% DEFINING GMC (Gaussian Mixture Copula) CLASS
classdef gmcdistribution
    
    properties (SetAccess = private)
        mu     % k x d mean vectors
        Sigma  % d x d x k  covariance matrices
        PComponents  % k x 1  mixing proportion        
    end
    
     properties (SetAccess = public)
        nLogL  % negative log likelihood
        numParams % number of parameters
        N      % Number of Data points used for fitting the model
    end
    
    properties (Dependent = true)
        BIC    % Bayesian Information Criterion
        NComponents % number of components
        NDimensions
    end
    
    methods
        % CONSTRUCTOR METHOD
        % Constructs a gmcdistribution object for a given parameter set.
        function obj = gmcdistribution(mu,Sigma,PComponents,varargin)  
            
            p = inputParser;
            p.addParamValue('num_dims',2);
            p.addParamValue('num_components',1);
            p.parse(varargin{:});
           
            if ~isempty(mu)
                obj.mu = mu;
                obj.Sigma = Sigma;
                obj.PComponents = PComponents;
            else
                obj.mu = zeros(1,p.Results.num_dims);
                obj.Sigma = zeros(p.Results.num_dims,p.Results.num_dims,p.Results.num_components);
                obj.PComponents = zeros(1,p.Results.num_components);
            end
        end
        
        % getting BIC
        function value = get.BIC(obj)
            value = 2*obj.nLogL + obj.numParams*log(obj.N);            
        end
        % getting number of components
        function value = get.NComponents(obj)
            value = length(obj.PComponents);            
        end
        % getting number of components
        function value = get.NDimensions(obj)
            value = size(obj.mu,2);            
        end
    end
    
end