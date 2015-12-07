%===================================================
function [indout] = ellipse(z, alpha)
%===================================================
% Drawing ellipse
% SYNOPSIS: ELLIPSE(X0, R, alpha)
% inputs:
%    X0 = Coordinates of the ellipse's center (2x1)
%    R  = A covariance positive (2x2) matrix
%    alpha  = confidence level 
% outputs:
% areaMC : area using monte-carlo
% area : area of the confidence ellipse at 100alpha%
%===================================================

z     = z(:);
N     = length(z);
meanz = nanmean(z);
zc    = z-ones(N,1)*meanz;
HH    = [real(zc) imag(zc)];
R     = HH'*HH/N;
rizc  = [real(zc), imag(zc)];

X0    = [real(meanz);imag(meanz)];

c     = -2*log(1-alpha);
Nc     = 100; 
theta = (0:Nc) * (2*pi) ./ Nc ;
Y     = sqrt(c)*[cos(theta);sin(theta)];
Fm1   = sqrtm(R);
X     = diag(X0)*ones(2,Nc+1)+Fm1*Y;
plot(X(1,:),X(2,:),'r','linewidth',1);
grid on
valp   = eig(R);
area   = sqrt(prod(valp))*c*pi;
hold on
plot(real(z),imag(z),'x')
hold off
indout = zeros(N,1);
for ii=1:N
    indout(ii) = rizc(ii,:) * Fm1 *rizc(ii,:)'<c;
end
% %=====
% G=100000;
% invR=inv(R);
% AA=sqrt(max(valp*c));
% points=2*AA*(rand(G,2)-1/2);
% ss=zeros(G,1);
% for ii=1:G
%     ppaux=points(ii,:);
%     ss(ii)=ppaux*invR*ppaux'<c;
% end
% areaMC = sum(ss)*4*AA*AA/G;
%=====================================================