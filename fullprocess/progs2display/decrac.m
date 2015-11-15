function locations = decrac(x)
%===================================
taillemod=500;
ordre=30; % ordre du modele
seuil=7;  % regle des 3 sigma
taillefen=taillemod+2*(1+ordre);
recouvrement=3*(1+ordre);
decalage=taillefen-recouvrement+1;
RACPIS2 = sqrt(pi/2);
%-----
x=[zeros(1,ordre+1) x']; % on met des "0" pour le filtrage initial
taillefichier = length(x);
nbfen=fix(1+(taillefichier-taillefen)/decalage);
locations = nan(nbfen,1);
%
for indfen=0:nbfen-1,
    iptlect=indfen*decalage;
    xsig=x(iptlect+1:iptlect+taillefen);
    %%%calculdetec
    Ri = zeros(ordre+1,1);
    for ip=1:ordre+1
        Ri(ip) = xsig(1:taillefen-ip+1)*xsig(ip:taillefen)'/taillefen;
    end
    AiL=levinson(Ri);
    res_int=filter(AiL,1,xsig);
    res_int=res_int(ordre+1:taillefen);
    FAi=AiL(ordre+1:-1:1);
    residuel=filter(FAi,1,res_int);
    residueltr=residuel(2*ordre:end);
    absresiduel=abs(residueltr);
    sigma=RACPIS2*std(residueltr);%(absresiduel*ones(taillemod,1))/taillemod;
    [valmax, posmax]=max(absresiduel);
    if (valmax>=sigma*seuil)
        locations(indfen+1)=iptlect+posmax;
    end
end