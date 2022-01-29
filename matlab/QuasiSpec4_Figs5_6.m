clear
% Calculate and plot quasiparticle propagator and fit for mu/Delta =4

% Colors
col1=[0,68,102]/255; % UA donkerblauw
col2=[85,170,51]/255; % FTEW groen
col3=[136,17,51]/255; % UA Donkerrood
col4=[221,153,17]/255; % UA Oranje
col5=[170,34,136]/255; % FTIW Magenta
col6=[0,102,170]/255; % FWET blauw
col7=[255,204,0]/255; % FLW geel

% Filenames
ptFile='pertSelfEnBCS4F.mat';
scFile='SCenergyBCS4F.mat';
qpFile='FortSelfEnx4.mat';

% Perturbative solution
[kV,EK,kThPT,x0,xiPT,epPT]=loadPertEnergy(ptFile);

% Mean-field energy
e0=sqrt((kV.^2-x0).^2+1);
kV0=linspace(0,3.5,1001);
xik0=kV0.^2-x0;
e00=sqrt(xik0.^2+1);

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

% Data points for cleaner contour plot, separated by sector
[uKplot1,Kplot1,Zplot1,GCHplot1,~,~,~,~] = loadQuasiProp('FortSelfEnx4_SP1.mat'); 
[uKplot2,Kplot2,Zplot2,GCHplot2,~,~,~,~] = loadQuasiProp('FortSelfEnx4_SP2.mat'); 

sl=2;
zlims=[0,4];
lvllst=linspace(zlims(1),zlims(2),201);
figure
subplot(2,1,1)% Energy spectrum
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
    plot(kV,real(EK),'--','Color',col1,'LineWidth',1.5);
    plot(kV0,e00,'k:','LineWidth',1.5);

    % Threshold energies
    plot(uK,th12,'-','Color',[.5,.5,.5],'LineWidth',1.5);
    plot(uK,th13,'-','Color',[.5,.5,.5],'LineWidth',1.5);

    % Fit solutions
    plot(uK,real(zSC),'-','Color',col2,'LineWidth',1.5);

    % fit error bars
    plot(uK,real(zSC)+sl*real(sigZ) ...
        ,'-','Color',1.2*col2,'HandleVisibility','off');
    plot(uK,real(zSC)-sl*real(sigZ) ...
        ,'-','Color',1.2*col2,'HandleVisibility','off');
    
    % % k_th-lines
    % plotKline(kTh,[0,6],col2);
    % plotKline(kThPT,[0,5],col1);

    % Limits
    xlim([0,3.5]); ylim([0.9,6]); caxis(zlims);
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
        '$z_k$', ...
    },'Interpreter','latex','FontSize',14,'Location','southwest')
hold off
% Damping rate
subplot(2,1,2)
hold on
    % Perturbative
    plot(kV,-2*imag(EK),'--','Color',col1,'LineWidth',1.5);

    % Fit results
    plot(uK,-2*imag(zSC),'-','Color',col2,'LineWidth',1.5);

    % % k-lines
    % plotKline(kTh,[0,6],col2);
    % plotKline(kThPT,[0,6],col1);

    % fit error bars
    plot(uK,-2*imag(zSC)+2*sl*imag(sigZ) ...
        ,'-','Color',1.2*col2,'HandleVisibility','off');
    plot(uK,-2*imag(zSC)-2*sl*imag(sigZ) ...
        ,'-','Color',1.2*col2,'HandleVisibility','off');
    
    % Limits
    xlim([0,13.6]); ylim([0,2]);

    % Labels
    colorbar
    xlabel('$\hbar k/\sqrt{2m\Delta}$','Interpreter','latex','FontSize',14);
    ylabel('$\hbar \Gamma/\Delta$','Interpreter','latex','FontSize',14);
hold off



% Bogoliubov factors Mean field
Uk0=sqrt((1+xik0./e00)/2);
Vk0=sqrt((1-xik0./e00)/2);

% Bogoliubov factors perturbative
UkPT=sqrt((1+xiPT./epPT)/2);
VkPT=sqrt((1-xiPT./epPT)/2);

% Bogoliubov factors selfconsistent
UkSC=zeros(size(uK));
VkSC=zeros(size(uK));

% Self-energy
PP=load(qpFile); % x0 K uK Z SPP SMM SPM th12 th13

for ik=1:length(uK)
    % Select values
    zV=Z(ik,:);
    SPP=PP.SPP(ik,:);
    SMM=PP.SMM(ik,:);
    SPM=PP.SPM(ik,:);
    % Remove NaN values
    SPP=SPP(~isnan(zV));
    SMM=SMM(~isnan(zV));
    SPM=SPM(~isnan(zV));
    zV=zV(~isnan(zV));

    % Select unique values
    [zV,ia,~]=unique(zV);

    % Interpolate to Re[z_k+i 0^+]
    spp=interp1(zV,SPP(ia),real(zSC(ik)));
    smm=interp1(zV,SMM(ia),real(zSC(ik)));
    spm=interp1(zV,SPM(ia),real(zSC(ik)));

    % Bogoliubov coefficients
    XICH=uK(ik).^2-x0+(smm-spp)/2;
    EPCH=sqrt(XICH.^2+(1-spm).^2);
    UkSC(ik)=sqrt((1+XICH./EPCH)/2);
    VkSC(ik)=sqrt((1-XICH./EPCH)/2);

end

% Plot
figure
subplot(2,1,1)
hold on
    plot(uK,real(UkSC),'LineWidth',1.5);
    plot(kV,real(UkPT),'--','LineWidth',1.5);
    plot(kV0,Uk0,'k:','LineWidth',1.5);
    xlim([0,3.5]);
    xlabel('$\hbar k/\sqrt{2m\Delta}$','Interpreter','latex');
    ylabel('$\mathcal{U}$','Interpreter','latex');
    legend({'$\mathrm{Re} \mathcal{U}(k,\mathrm{Re} z_{k})$', ...
        '$\mathrm{Re} \mathcal{U}(k,\mathrm{Re} \epsilon_k+i 0^+)$', ...
        '$U_k$'},'Interpreter','latex');
hold off
subplot(2,1,2)
hold on
    plot(uK,real(VkSC),'LineWidth',1.5);
    plot(kV,real(VkPT),'--','LineWidth',1.5);
    plot(kV0,Vk0,'k:','LineWidth',1.5);
    xlim([0,3.5]);
    xlabel('$\hbar k/\sqrt{2m\Delta}$','Interpreter','latex');
    ylabel('$\mathcal{V}$','Interpreter','latex');
    legend({'$\mathrm{Re} \mathcal{V}(k,\mathrm{Re} z_{k})$', ...
        '$\mathrm{Re} \mathcal{V}(k,\mathrm{Re} \epsilon_k+i 0^+)$', ...
        '$V_k$'},'Interpreter','latex');
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