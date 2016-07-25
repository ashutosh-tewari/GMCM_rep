% Author: Ashutosh Tewari (tewaria@utrc.utc.com)
% Affiliation: Decision Support & Machine Intelligence
%              United Technologies Research Center
%              East Hartford, CT: 06118
% Code based on the paper 'Parametric Characterization of MultiModal Densities with non-Gaussian
% Nodes', A. Tewari, M.J. Giering, A. Raghunathan, OEDM workshop, ICDM 2011 (paper included in the package) 

% DEFINING GMCM (Gaussian Mixture Copula Model) CLASS
classdef GMCM < handle
    
    properties (SetAccess = private)
        marginals   % cell array with marginal densities
        gmc  % Gaussian mixture copula (gmc) density
    end
    
    properties (SetAccess = public)
        model_setting % parameters indicating different choices
        nLogL  % negative log likelihood
        numParams % number of parameters
        NDataSamples      % Number of Data points used for fitting the model
        NComponents     % number of components to be used in gmc density
        NDimensions     % number of dimensons
    end
    
    properties (Dependent = true)
        BIC    % Bayesian Information Criterion
    end
    
    methods
        % CONSTRUCTOR METHOD
        % Constructs a gmcdistribution object for a given parameter set.
        function obj = GMCM(varargin)  
            
            p = inputParser;
            p.addParamValue('marginal_type','nonparam');
            p.addParamValue('useAnalyticalGradient',0);
            p.addParamValue('useParallel',0);
            p.addParamValue('num_EM_iterations',100);
            p.addParamValue('num_GD_iterations',0);
            p.parse(varargin{:});
                        
            obj.model_setting.marginal_type = p.Results.marginal_type; % Specify if you want to use emperical or non-parametric CDF of marginals.
            obj.model_setting.useAnalyticalGradient = p.Results.useAnalyticalGradient; % Specify if you want to compute gradient  analytically or by Finite Difference.
            obj.model_setting.useParallel = p.Results.useParallel; % Specify if you want to use parallel computing toolbox for gradient computation.
            obj.model_setting.num_EM_iterations = p.Results.num_EM_iterations; % Max number of EM iterations
            obj.model_setting.num_GD_iterations = p.Results.num_GD_iterations; % Max number of Gradient Descent iterations

        end
        
        % getting BIC
        function value = get.BIC(obj)
            value = 2*obj.nLogL + obj.numParams*log(obj.NDataSamples);            
        end

    end
    
end