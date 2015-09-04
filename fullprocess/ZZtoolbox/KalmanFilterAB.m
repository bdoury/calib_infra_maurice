%==========================================================
function [numB, denA, Xtt, Ptt, Pttm1] = ...
    KalmanFilterAB(Y, X, Lya, Lxb, Lc, RW, mu0)
%==========================================================
% Kalman's filter
%========
% Inputs:
% RV : state covariance
% RW : observation covariance
% Y  : observations
% mu0 : initial state mean
% R0 : initial state covariance
%
% Outputs:
%   Xfilter_k : filtered state
%   Covfilter_k : covariance
%   p(X_k|y_{1:k}) \sim Gauss(Xfilter_k, Covfilter_k)
%==========================================================
RW = RW*eye(Lc);
k0=max(Lxb+Lya);

dimX                = Lya+Lxb;
T                   = length(Y);
Xtt                 = zeros(dimX,T); %Acols*Ycols
Ptt                 = zeros(dimX,dimX,T); %Acols*Acols%Ycols
Ptt(:,:,1)          = 0;
Pttm1               = zeros(dimX,dimX,T);

Xttm1              = mu0;
Xtt(:,k0)        = Xttm1;
Pttm1(:,:,1)     = eye(Lxb+Lya);
Ptt(:,:,k0)      = eye(Lxb+Lya);
Id=eye(dimX);
for k=k0+1:T-Lya-Lxb-Lc,
    if Lya==0
        Bxk = toeplitz(X(k:k+Lc-1),X(k:-1:k-Lxb+1));
        Bk  = Bxk;
    elseif Lxb==0
        Byk = toeplitz(Y(k-1:k+Lc-2),Y(k-1:-1:k-Lya));
        Bk  = -Byk;
    else
        Byk = toeplitz(Y(k-1:k+Lc-2),Y(k-1:-1:k-Lya));
        Bxk = toeplitz(X(k:k+Lc-1),X(k:-1:k-Lxb+1));
        Bk  = [-Byk Bxk];
    end
    Xttm1              = Xtt(:,k-1);
    Rttm1              = Ptt(:,:,k-1);
    cov_inov           = (Bk * Rttm1 * Bk') + RW;
    inov_k             = Y(k:k+Lc-1)- Bk*Xttm1;
    Kn                 = (Rttm1*Bk') / cov_inov;
    Xtt(:,k)           = Xttm1 + Kn * inov_k;
    Pttm1(:,:,k)       = Rttm1;
    Ptt(:,:,k)         = (Id - Kn*Bk)*Rttm1;
end
Xtt=Xtt(:,k0+1:T-Lya-Lxb-Lc);
Ptt=Ptt(:,:,k0+1:T-Lya-Lxb-Lc);
Pttm1=Pttm1(:,:,k0+1:T-Lya-Lxb-Lc);
denA = [1 ; Xtt(2:Lya,end)];
numB = Xtt(Lya+1:Lya+Lxb,end);

%==========================================================
