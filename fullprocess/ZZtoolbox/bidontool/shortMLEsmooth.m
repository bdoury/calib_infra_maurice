function [x,Jval] = shortMLEsmooth(R,g,sigma2)
[x,Jval] = fminbnd(@(x) J(x,R,g,sigma2),0,100);
function valJ = J(x,R,g,sigma2)
S = [x*x x;x 1]+sigma2*[1 0;0 g];
valJ = real(log(det(S))+trace(S\R));
%=====================================================