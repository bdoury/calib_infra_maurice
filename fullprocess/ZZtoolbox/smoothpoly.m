function [xe,ye] = smoothpoly(x,y,P,degs,N)
% degs = (0:1.5:7);
x=x(:);
y=y(:);
xp=x(not(isnan(y)));
yp=y(not(isnan(y)));
yp0=yp(not(xp==0));
xp0=xp(not(xp==0));
xe = zeros(N,1);
ye = zeros(N,1);
logspace_xe = logspace(log10(xp0(1)),log10(xp0(end)),N)';

NonP = round(N/P);
for ip=0:P-1
    ide1=ip*NonP+1;
    ide2=min([ide1+NonP,N]);
    xe(ide1:ide2) = logspace_xe(ide1:ide2);

    ind1_ip = find(xp0>xe(ide1),1,'first');
    ind2_ip = find(xp0<xe(ide2),1,'last');
    xp_ip = xp0(ind1_ip:ind2_ip);
    yp_ip = yp0(ind1_ip:ind2_ip);
    Hp_ip = exp(log(xp_ip)*degs)+0.00000001*eye(length(xp_ip),length(degs));
    alpha_ip=Hp_ip\yp_ip;
   
    He_ip =  exp(log(xe(ide1:ide2))*degs);
    ye(ide1:ide2)=He_ip*alpha_ip;
end
ye = filtfilt(ones(3,1)/3,1,ye);