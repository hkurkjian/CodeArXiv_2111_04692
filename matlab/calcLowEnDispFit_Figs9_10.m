clear
% Fit the low-energy dispersion to extract effective gap, mass, and location of the minimum

flrd='SCenergy';
flrdPT='pertSelfEn';
% Values of mu/Delta
xVec=[0.8604366861256786,1.5,2.1,3,4,5,6,8,10,15,20,25,35,45,60,100,300,1000,5000];
xName=["Uni","1.5","2.1","3","4","5","6","8","10","15","20","25","35","45","60","100","300","1000","5000"];

%%%% Selfconsistent

% Initialize parameters
GapSC=zeros(size(xVec));
kmSC=zeros(size(xVec));
e0SC=zeros(size(xVec));
meffSC=zeros(size(xVec));

% Fit to minimum of the eigenenergy
for ix=1:length(xVec)

    % data file
    if strcmp(xName(ix),"Uni")
        lFile=[flrd,char(xName(ix))];
    else
        lFile=[flrd,'BCS',char(xName(ix))];
    end
    % Load data file
    if isfile([lFile,'F.mat'])
        F=load([lFile,'F.mat']);
    elseif isfile([lFile,'.mat'])
        F=load([lFile,'.mat']);
    else
        error(['No file found for x0=',char(xName(ix))]);
    end

    % See if zk is saved as zk or zV
    try
        F.zk;
    catch
        F.zk=F.zV;
    end

    % Remove NaN's
    kV=F.kV(~isnan(F.zk));
    zk=F.zk(~isnan(F.zk));

    % Make columns
    kV=kV(:);
    zk=zk(:);

    % Gap and location of minimum by interpolation
    fSCEn=@(k) interp1(kV,zk,k,'makima');
    [kmSC(ix),GapSC(ix)]=fminbnd(fSCEn,kV(1),kV(end));
    
    % Linear fit
    kCloseMin=(kmSC(ix)-.1./kmSC(ix)<=kV & kV<=kmSC(ix)+.1./kmSC(ix));
    t=((kV(kCloseMin)-kmSC(ix)).^2)';

    X=[ones(size(t));t]';
    Y=(zk(kCloseMin));
    B=X\Y;
    e0SC(ix)=B(1);
    meffSC(ix)=B(2).^(-1);

end

%%%% Perturbative

% Initialize parameters
GapPT=zeros(size(xVec));
kmPT=zeros(size(xVec));
e0PT=zeros(size(xVec));
meffPT=zeros(size(xVec));

% Fit to minimum of the eigenenergy
for ix=1:length(xVec)
    x0=xVec(ix);

    % data file
    if strcmp(xName(ix),"Uni")
        lFile=[flrdPT,char(xName(ix))];
    else
        lFile=[flrdPT,'BCS',char(xName(ix))];
    end
    % Load data file
    if isfile([lFile,'F.mat'])
        F=load([lFile,'F.mat']);

        [F.kV,iU,~]=unique(F.kV);

        Spp=F.rSpp(iU)+1i*F.iSpp(iU);
        Smm=F.rSmm(iU)+1i*F.iSmm(iU);
        Spm=F.rSpm(iU)+1i*F.iSpm(iU);
    elseif isfile([lFile,'.mat'])
        F=load([lFile,'.mat']);

        Spp=(F.Sig12{1,1}-1i*pi*F.iSig12{1,1} ...
            +F.Sig03{1,1} ...
            +F.Sig13{1,1}+1i*F.iSig13{1,1} ...
            +F.Sig04{1,1});
        Smm=(F.Sig12{2,2}-1i*pi*F.iSig12{2,2} ...
            +F.Sig03{2,2} ...
            +F.Sig13{2,2}+1i*F.iSig13{2,2} ...
            +F.Sig04{2,2});
        Spm=(F.Sig12{1,2}-1i*pi*F.iSig12{1,2} ...
            +F.Sig03{1,2} ...
            +F.Sig13{1,2}+1i*F.iSig13{1,2} ...
            +F.Sig04{1,2});
    else
        disp(['No perturbative file found for x0=',char(xName(ix))]);
        GapPT(ix)=NaN;
        kmPT(ix)=NaN;
        e0PT(ix)=NaN;
        meffPT(ix)=NaN;
        continue
    end

    % Perturbative energy
    F.zk=real(-Spp-Smm+sqrt((2*F.kV.^2-2*x0-Spp+Smm).^2+4*(1-Spm).^2))/2;

    % Remove NaN's
    kV=F.kV(~isnan(F.zk));
    zk=F.zk(~isnan(F.zk));

    % Make columns
    kV=kV(:);
    zk=zk(:);

    % Gap and location of minimum by interpolation
    fPTEn=@(k) interp1(kV,zk,k,'makima');
    [kmPT(ix),GapPT(ix)]=fminbnd(fPTEn,kV(1),kV(end));
    
    % Linear fit
    kCloseMin=(kmPT(ix)-.1./kmPT(ix)<=kV & kV<=kmPT(ix)+.1./kmPT(ix));
    t=((kV(kCloseMin)-kmPT(ix)).^2)';
    X=[ones(size(t));t]';
    Y=(zk(kCloseMin));
    B=X\Y;
    e0PT(ix)=B(1);
    meffPT(ix)=B(2).^(-1);

end

% Mean field
ainV0=linspace(-5.5,0,1001);
xVec0=xkF(ainV0,5);
aVec0=ainV0.^(-1);

e00=ones(size(xVec0));
km0=sqrt(xVec0);
meff0=xVec0.^(-1)/2;

% 1/kF a 
D0=Delta(xVec);
ainV=kFx(xVec);
aVec=ainV.^(-1);


figure('OuterPosition',[2851,-127,915,1070])
subplot(2,1,1)
    plot(ainV0,e00,'k:',ainV,e0SC,'o-',ainV,e0PT,'x--','LineWidth',1.5);
    xlim([-5.5,0]); ylim([0.8,1.07])
    xlabel('$1/k_\mathrm{F} a$','Interpreter','latex','FontSize',14);
    ylabel('$\epsilon^\ast/\Delta$','Interpreter','latex','FontSize',14);
subplot(2,1,2)
    plot(ainV0,meff0.*(2*xVec0),'k:',ainV,meffSC.*(2*xVec),'o-',ainV,meffPT.*(2*xVec),'x--','LineWidth',1.5);
    xlim([-5.5,0]);
    xlabel('$1/k_\mathrm{F} a$','Interpreter','latex','FontSize',14);
    ylabel('$m_\mathrm{eff}/m$','Interpreter','latex','FontSize',14);
    legend({'Mean-field','Self-consistent','Perturbative'},'Interpreter','latex','Location','NorthWest','FontSize',14)


% Quadratic shift
sh2=kmSC.^2-xVec;
sh2F=sh2.*D0;

% Theoretical
kfaHS=linspace(-2,0,101);
HarShift=-4/3/pi*kfaHS;
GalShift=-4/3/pi*kfaHS-4/15/pi^2*(11-2*log(2))*kfaHS.^2;

figure('OuterPosition',[1914,-6,974,809])
hold on
    plot(kfaHS,HarShift,'k:',kfaHS,GalShift,':',aVec,sh2F,'o-','LineWidth',1.5);
    xlabel('$k_\mathrm{F} a$','Interpreter','latex','FontSize',13);
    ylabel('$((k^\ast_m)^2-(k_m^{(0)})^2)/k_\mathrm{F}^2 $','Interpreter','latex','FontSize',13);
    xlim([-1.5,0]); ylim([0,.4]);
    legend({ ...
        'Hartree [$-\frac{4}{3\pi} k_\mathrm{F} a$]', ...
        'Galitskii[$-\frac{4}{3\pi} k_\mathrm{F} a - \frac{4(11-2\log 2)}{15\pi^2} (k_\mathrm{F} a)^2$]', ...
        'Selfconsistent', ...
        }, ...
    'Interpreter','latex','FontSize',12)
hold off