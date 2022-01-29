function [uK,K,Z,GCH,th12,th13,kTh,x0] = loadQuasiProp(filename)
% Calculate the quasiparticle propagator \check{G}
% from the values stored in filename (.mat file)
% output:
%   uK    -> unique values of hbar k/sqrt(2 m Delta)
%   K     -> matrix hbar k/sqrt(2 m Delta)
%   Z     -> matrix z_k/Delta
%   GCH   -> 2x2 cell with GCH{1,1}=G_{+,+}, GCH{1,2}=G_{+,-}, GCH{2,2}=G_{-,-}
%   th12  -> epsilon_th^(1->2)/Delta
%   th13  -> epsilon_th^(1->3)/Delta
%   kTh   -> Threshold values where zk=epsilon_th (1x4 vector)
%   x0    -> number mu/Delta

    % Load file
    F=load(filename); % x0 K uK Z SPP SMM SPM th12 th13
    K=F.K;
    uK=F.uK;
    Z=F.Z;
    x0=F.x0;
    th12=F.th12;
    th13=F.th13;

    % Calculate inverse Green's function
    invGrFun=cell(2,2);
    invGrFun{1,1}=-Z+K.^2-x0-F.SPP;
    invGrFun{2,2}=-Z-K.^2+x0-F.SMM;
    invGrFun{1,2}=1-F.SPM;

    % Generalized Bogoliubov matrix
    XICH=K.^2-x0+(F.SMM-F.SPP)/2;
    EPCH=sqrt(XICH.^2+(1-F.SPM).^2);
    UK=sqrt((1+XICH./EPCH)/2);
    VK=sqrt((1-XICH./EPCH)/2);

    % Quasiparticle basis
    GCHinv=cell(2,2);
    GCHinv{1,1}=(invGrFun{1,1}.*UK+invGrFun{1,2}.*VK).*conj(UK)+(invGrFun{1,2}.*UK+invGrFun{2,2}.*VK).*conj(VK);
    GCHinv{2,2}=(invGrFun{2,2}.*UK-invGrFun{1,2}.*VK).*conj(UK)-(invGrFun{1,2}.*UK-invGrFun{1,1}.*VK).*conj(VK);
    GCHinv{1,2}=(invGrFun{1,2}.*UK+invGrFun{2,2}.*VK).*conj(UK)-(invGrFun{1,1}.*UK+invGrFun{1,2}.*VK).*conj(VK);
    
    % Determinant
    DETGCH=GCHinv{1,1}.*GCHinv{2,2}-GCHinv{1,2}.^2;

    % Invert inverse Green's function
    GCH=cell(2,2);
    GCH{1,1}= GCHinv{2,2}./DETGCH;
    GCH{2,2}= GCHinv{1,1}./DETGCH;
    GCH{1,2}=-GCHinv{1,2}./DETGCH;

    % Threshold values
    if isfile([filename(1:end-4),'_kth.mat'])
        % Load from file
        T=load([filename(1:end-4),'_kth.mat']);
        kTh=zeros(1,4);
        kTh(1)=T.kTh13M;
        kTh(2)=T.kTh12M;
        kTh(3)=T.kTh12P;
        kTh(4)=T.kTh13P;
    else
        error(['Threshold values should be stored in ',filename(1:end-4),'_kth.mat']);
    end

end