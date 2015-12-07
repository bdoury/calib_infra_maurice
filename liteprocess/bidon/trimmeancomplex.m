function trimmedz = trimmeancomplex(z,apercent)

[ra,co] = size(z);
trimmedz = nan(ra,co);

for ira=1:ra 
    indout   = quadform(z(ira,:),apercent);
    trimmedz(ira,indout==1) = z(ira,indout==1);
end
%===================================================
function indout = quadform(z, apercent)
%===================================================
c     = -2*log(1-apercent);
z     = z(:);
N     = length(z);
meanz = nanmean(z);
zc    = z-ones(N,1)*meanz;
HH    = [real(z) imag(z)];
R     = nancov(HH);
rizc  = [real(zc), imag(zc)];
if or(sum(any(isnan(R)))>0,sum(any(isinf(R)))>0)
    indout = zeros(N,1);
elseif rank(R)==2
    Fm1   = sqrtm(R);
    % valp   = eig(R);
    % area   = sqrt(prod(valp))*c*pi;
    indout = zeros(N,1);
    for ii=1:N
        indout(ii) = rizc(ii,:) * Fm1 *rizc(ii,:)'<c;
    end
else
    indout = zeros(N,1);
end
%=====================================================