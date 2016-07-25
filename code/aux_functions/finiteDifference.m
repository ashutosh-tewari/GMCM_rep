function grad = finiteDifference(fun,x0)

delta = 1E-8;

grad = nan(numel(x0),1);
for i = 1:numel(x0)
    x_plus = x0;x_minus = x0;
    x_plus(i) = x_plus(i) + delta;
    x_minus(i) = x_minus(i) - delta;
    grad(i) = (fun(x_plus) - fun(x_minus))/(2*delta);    
end