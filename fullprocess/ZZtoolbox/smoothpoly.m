function [xe,ye] = smoothpoly(x,y,degs,N)
% degs = (0:1.5:7);
x=x(:);
y=y(:);
xp=x(not(isnan(y)));
yp=y(not(isnan(y)));
xp=xp(2:end);
yp=yp(2:end);
Hp = exp(log(xp)*degs)+00.0001*eye(length(xp),length(degs));
alpha=Hp\yp;
% xe=linspace(xp(1),xp(end-1),N)';
xe = logspace(log10(xp(1)),log10(xp(end)),N)';

He =  exp(log(xe)*degs);
ye=He*alpha;
