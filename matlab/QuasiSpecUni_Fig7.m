clear 
% Calculate and plot quasiparticle propagator and fit at unitarity

% Colors
col1=[0,68,102]/255; % UA donkerblauw
col2=[85,170,51]/255; % FTEW groen
col3=[136,17,51]/255; % UA Donkerrood
col4=[221,153,17]/255; % UA Oranje
col5=[170,34,136]/255; % FTIW Magenta
col6=[0,102,170]/255; % FWET blauw
col7=[255,204,0]/255; % FLW geel

% Perturbative solution
[kV,EK,kThPT,x0]=loadPertEnergy('pertSelfEnUniF.mat');
e0=sqrt((kV.^2-x0).^2+1);

% Selfconsistent solution
SC=load('SCenergyUniF.mat'); % kV zV th12 th13

% Quasiparticle propagator
[uK,K,Z,GCH,th12,th13,kTh,x0] = loadQuasiProp('FortSelfEnUni.mat'); 

% Fit Lorentzian in Different regions
aIn= 0.8805 - 0.0670i;
zIn= 1.7600 - 0.0071i;

% before kTh(3)=1.8157, only 1 pole
pt1=find(uK<kTh(3),1,'last');
uK1=uK(1:pt1);
[zSC1,aSC1,sigA1,sigZ1,rms1]=calcSpectrumFitUni( ...
    aIn,zIn, ...
    uK1,Z(1:pt1,:),GCH{1,1}(1:pt1,:), ...
    kTh,SC,th12(1:pt1),th13(1:pt1), ...
    'Plot',[],'FiltRe',0.001,'Filt',1/2);

% Below th13, to k=3.05
kEND2=2.8;
aIn=[0.5425 - 0.2467i,-0.0511 - 0.1061i];
zIn=[1.8183 - 0.0129i,1.8252 - 0.0008i];
kTh(4)=5;
pt2=find(uK==kEND2);
uK2=uK(pt1+1:pt2);
dblPts=true(size(uK2));
[zSC2,aSC2,sigA2,sigZ2,rms2]=calcSpectrumFitUni( ...
    aIn,zIn, ...
    uK2,Z(pt1+1:pt2,:),GCH{1,1}(pt1+1:pt2,:), ...
    kTh,SC,th12(pt1+1:pt2),th13(pt1+1:pt2), ...
    'Plot',[],'Filt',1/2,'DblLz',dblPts,'FiltRe',0.1);

% Above th13, from k=1.95
kBEG3=1.95;
aIn=1.0015 - 0.3910i;
zIn=2.5302 - 0.4231i;
kTh(4)=1;
pt3=find(uK==kBEG3);
uK3=uK(pt3:end);
dblPts=true(size(uK3));
rgm=NaN*zeros(size(uK3));
rgm(uK3>2.5)=0;
[zSC3,aSC3,sigA3,sigZ3,rms3]=calcSpectrumFitUni( ...
    aIn,zIn, ...
    uK3,Z(pt3:end,:),GCH{1,1}(pt3:end,:), ...
    kTh,SC,th12(pt3:end),th13(pt3:end), ...
    'Plot',[],'Filt',1/2,'DblLz',dblPts,'Rgm',rgm);

myK=2.16;
ik0=find(uK==myK);
ik2=find(uK2==myK);
ik3=find(uK3==myK);

% Data points for cleaner contour plot, separated by sector
[uKplot1,Kplot1,Zplot1,GCHplot1,~,~,~,~] = loadQuasiProp('./PertEnSave/FortSelfEnUni_SP1.mat'); 
[uKplot2,Kplot2,Zplot2,GCHplot2,~,~,~,~] = loadQuasiProp('./PertEnSave/FortSelfEnUni_SP2.mat'); 

sl=2;
zlims=[0,4];
lvllst=linspace(zlims(1),zlims(2),201);
figure
subplot(2,1,1) % Energy spectrum
hold on
    % Spectral density
    [~,cD1]=contourf(Kplot1,Zplot1,imag(GCHplot1{1,1}));
    cD1.LevelList=lvllst;
    cD1.LineStyle='none';
    [~,cD2]=contourf(Kplot2,Zplot2,imag(GCHplot2{1,1}));
    cD2.LevelList=lvllst;
    cD2.LineStyle='none';
    cD2.HandleVisibility='off';
    
    % Mean-field & perturbative
    plot(kV,real(EK),'--','Color',col6,'LineWidth',1.5);
    plot(kV,e0,'k:','LineWidth',1.5);

    % Threshold energies
    plot(uK,th12,'-','Color',[.5,.5,.5],'LineWidth',1.5);
    plot(uK,th13,'-','Color',[.5,.5,.5],'LineWidth',1.5);

    % Fit solutions
    uKN=[uK1;uK2];
    eSC1=real([zSC1(:,1);zSC2(:,1)]);
    eSC3=real(zSC3(:,1));
    dSC1=-2*imag([zSC1(:,1);zSC2(:,1)]);
    dSC3=-2*imag(zSC3(:,1));

    plot(uKN(eSC1<3),eSC1(eSC1<3),'-','Color',col2,'LineWidth',1.5);
    plot(uKN(eSC1>3),eSC1(eSC1>3),':','Color',col2,'LineWidth',1.5,'HandleVisibility','off');
    plot(uK3(eSC3<3),eSC3(eSC3<3),':','Color',col3,'LineWidth',1.5,'HandleVisibility','off');
    plot(uK3(eSC3>3),eSC3(eSC3>3),'-','Color',col3,'LineWidth',1.5);
    plot(uK2,real(zSC2(:,2)),'-','Color',col6,'LineWidth',1.5);
    plot(uK3(~isnan(zSC3(:,1))),real(zSC3(~isnan(zSC3(:,1)),2)),'-','Color',col5,'LineWidth',1.5);

    % Limits
    xlim([0,3]); ylim([0.0,8.0]); caxis(zlims);
    colorbar; colormap(flipud(hot))

    % Labels
    xlabel('$\hbar k/\sqrt{2m\Delta}$','Interpreter','latex','FontSize',14);
    ylabel('$\varepsilon/\Delta$','Interpreter','latex','FontSize',14);

    % Legend
    legend({ ...
        '$\mathrm{Im}[\check{G}_{+,+}] \Delta$', ...
        '$z_k^{(2)}$', ...
        '$\epsilon_k$', ...
        '$\epsilon_\mathrm{th}^{(1\rightarrow 2)}$', ...
        '$\epsilon_\mathrm{th}^{(1\rightarrow 3)}$', ...
        '$z_k$ ($\varepsilon < \epsilon_\mathrm{th}^{1\rightarrow 3}$)', ...
        '$z_k$ ($\varepsilon > \epsilon_\mathrm{th}^{1\rightarrow 3}$)', ...
    },'Interpreter','latex','FontSize',14,'Location','northwest')
hold off
% Damping rate
subplot(2,1,2)
hold on
    % Perturbative
    plot(kV,-2*imag(EK),'--','Color',col6,'LineWidth',1.5);

    % Fit results
    plot(uKN(eSC1<3),dSC1(eSC1<3),'-','Color',col2,'LineWidth',1.5);
    plot(uKN(eSC1>3),dSC1(eSC1>3),':','Color',col2,'LineWidth',1.5);
    plot(uK3(eSC3<3),dSC3(eSC3<3),':','Color',col3,'LineWidth',1.5);
    plot(uK3(eSC3>3),dSC3(eSC3>3),'-','Color',col3,'LineWidth',1.5);
    
    % Limits
    xlim([0,3]); ylim([0,3.0]);

    % Labels
    % colorbar
    xlabel('$\hbar k/\sqrt{2m\Delta}$','Interpreter','latex','FontSize',14);
    ylabel('$\hbar \Gamma/\Delta$','Interpreter','latex','FontSize',14);
hold off



function plotKline(kTh,yLims,varargin)
    % Plot thin vertical lines

    if nargin<3
        col=[.5,.5,.5];
    elseif nargin<4
        col=varargin{1};
    end
    
    for ii=1:length(kTh)
        plot([kTh(ii),kTh(ii)],yLims,'-', ...
            'Color',col, ...
            'LineWidth',.5, ...
            'HandleVisibility','off');
    end
    
end