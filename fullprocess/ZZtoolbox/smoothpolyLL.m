function [xo,yo] = ...
    smoothpolyLL(xi,yi,P,degs,logtrain,logfit)
%==========================================================
% Smoothing with polynomials:
% Learning could be done on a part of the input values
%     Nan are removed
% Smoothing fit could be done on a given grid 
%     of values.
% Synopsis
%    [xo,yo] = smoothpolyLL(xi,yi,P,degs,logtrain,logfit)
% Rk: if FREQ_FLAG is set to 1, the half 
%    of the assumed frequency range is taken into account
% Possible improvement: adjust derivative at the junctions
%==========================================================
% Inputs:
%     (xi,yi): two arrays 1D
%     P : number of intervals with different fitting
%         each segment has the same number of values
%     degs : sequence of powers, ex: degs = (0:1.5:7)
%            NON INTEGERS CAN BE USED
%     logtrain: structure with 
%              .flag, 0 for linear frequency scaling 
%                     1 for log frequency scaling
%              .N number of input values for
%               training 
%     logfit: structure with 
%              .flag, 0 for linear frequency scaling 
%                     1 for log frequency scaling
%              .N number of output values 
% Ex:
%  [xo,yo] = smoothpolyLL(xi,yi,2,(0:1:6),logtrain,logfit)
%  with logtrain.flag = 0, logtrain.N = 300, 
%       logfit.flag = 1, logfit.N = 100
%  300 values linearly distruted on xi are divided 
%       in 2 segments for training
%  Polynomial order sequence are 0,1, ...,6
%  100 values log distruted for xo are used for 
%       performing yo
%==========================================================

FREQ_FLAG = 1;
% regularization factor
lambda  = 0.00000001;
%======
xi = xi(:);
yi = yi(:);
N1=length(xi);
if FREQ_FLAG
    xi=xi(1:fix(N1/2));
    yi=yi(1:fix(N1/2));
end
yinan = yi(not(isnan(yi)));
xinan = xi(not(isnan(yi)));

if logtrain.N>length(xinan)
    error('********* not enough data on input values ********')
end
if logtrain.flag
    [xt,yt] = logextract(xinan,yinan, logtrain.N);
else
    xt = xinan(1:logtrain.N);yt=yinan(1:logtrain.N);
end
xt=xt(:);
yt=yt(:);
xp=xt(not(isnan(yt)));
yp=yt(not(isnan(yt)));
yp0=yp(not(xp==0));
xp0=xp(not(xp==0));

N = logfit.N;
xo = zeros(N,1);
yo = zeros(N,1);
if logfit.flag
    space_xe = logspace(log10(xp0(1)),...
        log10(xp0(end)),N)';
else
    space_xe = linspace(xp0(1),xp0(end),N)';
end
NonP = round(N/P);
for ip=0:P-1
    ide1=ip*NonP+1;
    ide2=min([ide1+NonP,N]);
    xo(ide1:ide2) = space_xe(ide1:ide2);

    ind1_ip = find(xp0>xo(ide1),1,'first');
    ind2_ip = find(xp0<xo(ide2),1,'last');
    xp_ip = xp0(ind1_ip:ind2_ip);
    yp_ip = yp0(ind1_ip:ind2_ip);
    Hp_ip = exp(log(xp_ip)*degs)+...
        lambda*eye(length(xp_ip),length(degs));
    alpha_ip=Hp_ip\yp_ip;
   
    He_ip =  exp(log(xo(ide1:ide2))*degs);
    yo(ide1:ide2)=He_ip*alpha_ip;
end
yo = filtfilt(ones(3,1)/3,1,yo);
%=====================================================
function [Fo,Ho] = logextract(Fi,Hi, Nop)
% extract log scaling values wrt Fi

N=length(Hi);
Flog = logspace(log10(Fi(2)),log10(Fi(N)),Nop)';
Fo = zeros(Nop,1);
Ho = zeros(Nop,1);
Fo(1)=Fi(2);
Ho(1)=Hi(2);
ic1=2;
for in=2:Nop
    ic = find(Fi(ic1:N)>Flog(in),1,'first');
    if not(isempty(Fi(ic1+ic-1)))
        Fo(in)=Fi(ic1+ic-1);
        Ho(in)=Hi(ic1+ic-1);
        ic1=ic1+ic;
        inlast=in;
    end
end
Fo=Fo(1:inlast);
Ho=Ho(1:inlast);
[Fo, indFo] = sort(Fo);
Ho = Ho(indFo);
%==========================================================

