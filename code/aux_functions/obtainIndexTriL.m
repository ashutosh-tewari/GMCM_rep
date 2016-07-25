function idx = obtainIndexTriL(d)

idx = [];
for i = 1:d
    idx = [idx [(i-1)*d+i:i*d]];
end
