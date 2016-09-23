% This function plots the contours of likelihood values on the scatter plot
% of a 2 dimensional data.
function [xgrid,ygrid,Z] = plotFit(obj,data,varargin) %logLikelihoodVal,numMeshPoints,x_dim,y_dim)

    %INPUT: data (MxN, M=number of points, N=Dimension)
    %     : plo = binary variable (1 plot contour plot, 0 do not plot)
    %OUTPUT: xgrid,ygrid,Z ( Z contains the likelihood values of the points defined by xgrid and ygrid)

    %By Ashutosh Tewari (tewaria@utrc.utc.com)

    p = inputParser;
    p.addParamValue('plot_type','contour_plot');
    p.parse(varargin{:});
    
    for i = 1:length(obj.marginals)
        lower_limits(i) = obj.marginals{i}.icdf(0.00001);
        upper_limits(i) = obj.marginals{i}.icdf(0.99999);
    end
        
       


    figure;
    d = obj.NDimensions;
    subplot1(d,d,'Gap',[0.02 0.02],'XTickL','All','YTickL','All','FontS',10);

    % Scatter plot of the data overlayed with the contours of fitted probability density values
    for i = 1:d
        for j = 1:d
            subplot_idx = sub2ind([d,d],i,j);
            subplot1(subplot_idx);   hold on;
            x_limit = [lower_limits(i) upper_limits(i)];
            y_limit = [lower_limits(j) upper_limits(j)];
            if i == j                
                marginalPlots(data(:,i),obj.marginals{i});
                xlim(x_limit);
                continue;
            end
            bivar_data = data(:,[i j]);
            marginals = obj.marginals([i j]);
            marg_gmc = obj.gmc.obtainMarginal([i j]);
            bivariatePlots(bivar_data,marginals,marg_gmc,p.Results.plot_type,256);
            xlim(x_limit);
            ylim(y_limit);
        end
    end

end

%%
function marginalPlots(data,marginal_density)

    x = [marginal_density.icdf(0.001): marginal_density.BandWidth: marginal_density.icdf(0.999)];
    y = marginal_density.pdf(x);
    y = y/sum(y);
    counts = histc(data,x);
    counts = counts/sum(counts);
    bar(x,counts,'c','EdgeColor','b','linewidth',0.25);
    plot(x,y,'k','linewidth',2);

end


%%
function bivariatePlots(data,marg_densities,marg_gmc,plot_type,bin_num)

    if nargin <5
        bin_num = 512;
    end
   
    for j = 1:2    
        l_limit(j) = marg_densities{j}.icdf(0.001);
        u_limit(j) = marg_densities{j}.icdf(0.999);
        xmesh_inverse_space{j} = (l_limit(j):(u_limit(j)-l_limit(j))/(bin_num-1):u_limit(j));
        marginalLogLikelihood_grid{j} = log(marg_densities{j}.pdf(xmesh_inverse_space{j}));
        marginalCDFValues_grid{j} = marg_densities{j}.cdf(xmesh_inverse_space{j});
    end

    % creating  x,y mesh for plotting results
    [xgrid,ygrid] = meshgrid(xmesh_inverse_space{1},xmesh_inverse_space{2});

    % Getting the marginal log densitiy values
    [marg1,marg2] = meshgrid(marginalLogLikelihood_grid{1},marginalLogLikelihood_grid{2});

    % preparing the input matrix for copula density
    [xg,yg] = meshgrid(marginalCDFValues_grid{1},marginalCDFValues_grid{2});
    inputMatrix = [reshape(xg,numel(xg),1) reshape(yg,numel(yg),1)];
    clear xg yg;

    % getting the log copula density values
    copulaLogLikelihoodVals =  log(marg_gmc.pdf(inputMatrix)); %gmmCopulaPDF(inputMatrix,gmmObject,xmesh_inverse_space);
    Z = reshape(copulaLogLikelihoodVals,size(marg1,1),size(marg2,2));
    % combining the marginals with copula log densities
    Z = Z+marg1+marg2;

    % Getting the likelihood value from the log-likelihood.
    Z = exp(Z); 

    switch plot_type
        case 'contour_plot'
            contour(xgrid,ygrid,Z,60);
            plot(data(:,1),data(:,2),'ko','MarkerSize',2); % overlaying the data points on the plots
        case 'image_plot'
            imagesc([min(min(xgrid)) max(max(xgrid))],[min(min(ygrid)) max(max(ygrid))],Z);
            plot(data(:,1),data(:,2),'w.','MarkerSize',2); % overlaying the data points on the plots
    end  
    axis tight;

end
