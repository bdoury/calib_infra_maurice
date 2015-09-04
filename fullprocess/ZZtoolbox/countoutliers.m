function countoutliers(data,alpha)
%!===================================================!
%! SYNOPSIS: CONFIDENCEELLIPSE(X0, E, c)             !
%! Ellipse equation:                                 !
%!    (X-X0)'*(C^-1)*(X-X0)= c(alpha)                !
%!    X0 = coordinates of the ellipse's center (2x1) !
%!    C  = positive (2x2) matrix                     !
%!    alpha = confidence level in(0,1)               !
%!===================================================!
dataRI = [real(data) imag(data)];
mdata=mean(dataRI,1);
[N,bid]=size(data);
dataRIc=dataRI-ones(N,1)*mdata;
E=dataRIc'*dataRIc/N;

theta = 2*pi*(0:100)/100;
c = -2*log(1-alpha);
Y = sqrt(c)*[cos(theta);sin(theta)];
X = sqrtm(E)*Y;
plot(X(1,:)+X0(1),X(2,:)+X0(2),'g','linew',2);
