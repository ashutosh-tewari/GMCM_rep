function [u] = learn_marginals(obj,givenData,refData,varargin) 

[N,d] = size(givenData);

% parsing the input argument
p = inputParser;
p.addParamValue('type','nonparam');
p.addParamValue('numMeshPointVector',256);
p.parse(varargin{:});

if numel(p.Results.numMeshPointVector) == 1;
    numMeshPointVector = ones(d,1)*p.Results.numMeshPointVector;
end

u = nan(N,d);

switch p.Results.type    
    case 'nonparam'
        if isempty(refData)
            refData = givenData;
        end
        for i = 1:d
            temp_vec = refData(:,i);
            numMeshPoints = numMeshPointVector(i);
            [bandwidth,~,~]=kde(temp_vec,numMeshPoints); %obtained an optimal bandwidth from a separate algorithm
            obj.marginals{i} = fitdist(temp_vec,'kernel','width',bandwidth); % fitting matlab kernel density estimator using the bandwidth 
            u(:,i) = obj.marginals{i}.cdf(temp_vec);
        end

    case 'empirical'
    
        if isempty(refData)  % if reference data does not exisit.    
            N_refData = 0;
            data_combined = givenData;               
        else    
            N_refData = size(refData,1);
            data_combined = [givenData;refData];
        end
        givenData_idx = 1:N; 
        u_combined = zeros(N+N_refData,d); % initializing cdf value matrix.
        for j = 1:d
            [~,u_combined(:,j)] = ismember(data_combined(:,j),sort(data_combined(:,j)));
        end
        u_combined = u_combined/(N +  N_refData + 1);
        u = u_combined(givenData_idx,:);
    
    case 'param'
        prob_dists = {'gev','normal','tlocationscale'};
        for i = 1:d
            temp_vec = givenData(:,i);
            for j = 1:length(prob_dists)
               PD{j} = fitdist(temp_vec,prob_dists{j}); 
               BIC(j) = 2*PD{j}.NLogL + PD{j}.NumParams*log(length(temp_vec));
            end
            [~,min_idx]=min(BIC);
            obj.marginals{i} = PD{min_idx}; 
            u(:,i) = obj.marginals{i}.cdf(temp_vec);
        end
end


