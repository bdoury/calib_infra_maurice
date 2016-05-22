function m = multivarmedian(x)

[N,d] = size(x);
m = fminsearch(@(v) f(x,v), zeros(d,1));



function J = f(x,v)
[N,d]=size(x);
xminusv = x - ones(N,1)*v';
J = sum(sum(abs(xminusv)));