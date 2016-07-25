function fit(obj,data,K)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here

    [obj.NDataSamples,obj.NDimensions]= size(data);
    obj.NComponents = K;

    % learning marginal densities first
    [u_mat] = obj.learn_marginals(data,[],'type',obj.model_setting.marginal_type);

    % learning copula density then
    gmcObject = gmcdistribution([],[],[],'num_dims',obj.NDimensions,'num_components',K);
    obj.gmc = gmcObject.fit(u_mat,K,obj.NDimensions,obj.NDataSamples,'Start','randSample','algorithm','active-set','iteration_EM',obj.model_setting.num_EM_iterations,'iteration_GD',obj.model_setting.num_EM_iterations,'EM_rep',1);

end

