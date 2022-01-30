function [kV,EK,kTh,x0,xich,epch] = loadPertEnergy(filename)
% Load perturbative energy solution from file
% Also threshold values are calculated

    % Load file
    F=load(filename);
    kV=F.kV;
    x0=F.x0;

    % Self-energy
    [kV,iU,~]=unique(kV);

    Spp=F.rSpp(iU)+1i*F.iSpp(iU);
    Smm=F.rSmm(iU)+1i*F.iSmm(iU);
    Spm=F.rSpm(iU)+1i*F.iSpm(iU);

    % Perturbative energy
    EK=(-Spp-Smm+sqrt((2*kV.^2-2*x0-Spp+Smm).^2+4*(1-Spm).^2))/2;

    % Energy functionals
    xich=kV.^2-x0+(Smm-Spp)/2;
    epch=sqrt(xich.^2+(1-Spm).^2);

    % Threshold values
    kTh=zeros(1,4);
    [kTh(2),kTh(3),~]=kIntValsPole(x0);

    if x0>=0 && x0<=sqrt(8)
        kTh(1)=NaN;
        kTh(4)=sqrt(x0+sqrt(8));
    elseif x0>=0 
        kTh(1)=sqrt(x0-sqrt(8));
        kTh(4)=sqrt(x0+sqrt(8));
    else
        kTh(1)=NaN;
        kTh(4)=sqrt(x0+sqrt(9*x0^2+8));
    end

end