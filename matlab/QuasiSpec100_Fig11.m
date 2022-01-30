clear
% Calculate and plot quasiparticle propagator and fit for mu/Delta = 100

% Colors
col1=[0,68,102]/255; % UA donkerblauw
col2=[85,170,51]/255; % FTEW groen
col3=[136,17,51]/255; % UA Donkerrood
col4=[221,153,17]/255; % UA Oranje
col5=[170,34,136]/255; % FTIW Magenta
col6=[0,102,170]/255; % FWET blauw
col7=[255,204,0]/255; % FLW geel

% Filenames
ptFile='pertSelfEnBCS100F.mat';
scFile='SCenergyBCS100F_pasres.mat';
qpFile='FortSelfEnx100.mat';

% Perturbative solution
[kV,EK,kThPT,x0]=loadPertEnergy(ptFile);

% Mean-field energy
e0=sqrt((kV.^2-x0).^2+1);
kV0=linspace(0,14,1001);
e00=sqrt((kV0.^2-x0).^2+1);

% Selfconsistent solution
SC=load(scFile); % kV zV th12 th13

% Quasiparticle propagator
[uK,K,Z,GCH,th12,th13,kTh,~] = loadQuasiProp(qpFile); 

% Fit the quasiparticle spectrum to a Lorentzian
aIn=-1.0335 + 0.1221i;
zIn= 4.9767 - 0.3352i;
[zSC,aSC,sigA,sigZ,rms]=calcSpectrumFit( ...
    aIn,zIn, ...
    uK,Z,GCH{1,1}, ...
    kTh,SC,th12,th13, ...
    'Plot',[],'FiltRe',0.01,'Filt',1/2);

% Data points for cleaner contour plot
[uKplot,Kplot,Zplot,GCHplot,~,~,~,~] = loadQuasiProp('FortSelfEnx100_SP.mat'); 

sl=2;
zlims=[0,1];
xLi=[0,40];
lvllst=linspace(zlims(1),zlims(2),201);

kmPT=10.576482219763516^2;
kmSC=10.572795431047297^2;
km0=sqrt(x0)^2;

figure
subplot(2,1,1) % Energy spectrum
hold on
    % Spectral density
    [~,cD]=contourf(Kplot.^2-kmSC,Zplot,imag(GCHplot{1,1}));
    cD.LevelList=lvllst;
    cD.LineStyle='none';
    
    % Mean-field & perturbative
    plot(kV.^2-kmPT,real(EK),'--','Color',col1,'LineWidth',1.5);
    plot(kV0.^2-km0,e00,'k:','LineWidth',1.5);

    % Fit solutions
    plot(uK.^2-kmSC,real(zSC),'-','Color',col2,'LineWidth',1.5);
    
    % k-lines
    plotKline(kThPT.^2-km0,[0,6],col1);
    plotKline(kTh.^2-kmSC,[0,6],col2);

    % Limits
    xlim([-6,6]); ylim([0,6]); 
    caxis(zlims);
    colorbar; colormap(flipud(hot))

    % Labels
    xlabel('$\xi_k/\Delta$','Interpreter','latex','FontSize',14);
    ylabel('$\varepsilon/\Delta$','Interpreter','latex','FontSize',14);

    % Legend
    legend({ ...
        '$\mathrm{Im}[\check{G}_{+,+}] \Delta$', ...
        '$z_k^{(2)}$', ...
        '$\epsilon_k$', ...
        '$z_k$', ...
    },'Interpreter','latex','FontSize',14,'Location','southwest')
hold off
% Damping rate
subplot(2,1,2)
hold on
    % Perturbative
    plot(kV.^2-km0,-2*imag(EK),'--','Color',col1,'LineWidth',1.5);

    % Fit results
    zSC(229)=real(zSC(229))+1i*imag(zSC(230))/2;
    zSC(imag(zSC)>0)=real(zSC(imag(zSC)>0));
    plot(uK.^2-kmSC,-2*imag(zSC),'-','Color',col2,'LineWidth',1.5);

    % k-lines
    plotKline(kThPT.^2-km0,[0,6],col1);
    plotKline(kTh.^2-kmSC,[0,6],col2);
    
    % Limits
    xlim([-6,6]); ylim([0,0.015]);

    % Labels
    colorbar
    xlabel('$\xi_k/\Delta$','Interpreter','latex','FontSize',14);
    ylabel('$\hbar \Gamma/\Delta$','Interpreter','latex','FontSize',14);
hold off


function plotKline(kTh,yLims,varargin)
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
