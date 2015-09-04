function sigout = fbank(sigin,filtercharact,Fs_Hz)



%==== structure: Px1
% filtercharact.Norder
% filtercharact.Wlow
% filtercharact.Whigh
%===========================


P = length(filtercharact);
[T,d]     = size(sigin);
sigout    = zeros(T,d,P);
for ifilter = 1:P
    [filnum,filden] = butter(filtercharact(ifilter).Norder,...
        2*[filtercharact(ifilter).Wlow filtercharact(ifilter).Whigh]/Fs_Hz);
    sigout(T,d,ifilter) = filter(filnum,filden,signin);
end 