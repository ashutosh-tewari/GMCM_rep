function objs = obtainGMMFit(data)

numComp_max = 10;
opt = statset('MaxIter',200);

for i = 1:numComp_max
    objs{i} = gmdistribution.fit(data,i,'replicate',10,'regularize',10E-4,'Options',opt);
    BIC(i) = objs{i}.BIC
end

