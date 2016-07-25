clear;clc;close all

% Adding paths.
addpath(genpath(pwd));

% generating synthetic data
numSamples=1000;numDim=4;numComp=5;plo=1;
example_data = genData(numSamples,numDim,numComp,plo);

%%% loading presaved data
% load synthetic_data_UAI;
% example_data=data;

% standardizing the data
data = zscore(example_data);

% instantiating the model
gmcmObj = GMCM('marginal_type','nonparam','num_EM_iterations',100);

% fitting the model to data
num_components = 5;
gmcmObj.fit(data,num_components);

% plotting results
gmcmObj.plotFit(data,'plot_type','contour_plot');

