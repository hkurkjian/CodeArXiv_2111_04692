function [zSC,aSC,sigA,sigZ,rms]=fitLorentzQSpec(aIn,zIn,uK,Z,QG,kTh,SC,th12,th13,varargin)
% fit the quasiparticle Green’s function to a Lorentzian resonance
%    G(k,epsilon + i 0^+) = a/(zk - epsilon)
% 
% Input variables:
%   aIn  -> Starting value for residue
%   zIn  -> Starting value for eigenenergy
%   uK   -> Vector with values of k
%   Z    -> Matrix with values of epsilon
%   QG   -> Matrix with values of the Quasiparticle Green's function
%   kTh  -> Vector with threshold values of k, where z_{k_th}=e_th(k_th)
%   SC   -> Structure containing the self-consistent eigenenergy solution zk
%           below e_th^{1->2}, and its k-values
%   th12 -> Vector with values of epsilon_th^{1->2}
%   th13 -> Vector with values of epsilon_th^{1->3}
%   
% Optional name-value parameters:
%   'Plot', [k1,k2]  ->  Plot the quasiparticle Green's function and its fit
%                        in function of the energy for values k1 < k < k2
%                        Default option does not plot anything
%   'Filt', fltr     ->  Filter data points for fit to only contain values
%                        up to fltr * max(QG)
%                        Default is set to fltr=1/2
%   'FiltRe', fR     ->  Filter data points for fit below epsilon_th^{1->2}
%                        by selecting values between zk-fR < epsilon < zk+fR
%                        with zk the self-consistent eigenenergy solution
%                        Default is set to fR=0.05
%   'SepReg', SR     ->  Boolean indicating whether values of the quasiparticle
%                        Green's function be split according to the different
%                        sectors defined in Fig. 4
%                        Default is set to SR=true
% 
% Output variables:
%   zSC   ->  Vector containing fitted values of the eigenenergy
%   aSC   ->  Vector containing fitted values of the residues
%   sigA  ->  Vector containing the standard deviation of the residues 
%   sigZ  ->  Vector containing the standard deviation of the eigenenergy
%   rms   ->  Vector containing the root-mean-square error

    % Extra inputs
    p=inputParser; % Generate input parser
    addParameter(p,'Plot',[]); % Optional parameter: plot region
    addParameter(p,'Filt',1/2); % Optional parameter: filter for fit
    addParameter(p,'FiltRe',0.05); % Optional parameter: filter for fit below th12
    addParameter(p,'SepReg',true); % Optional parameter: separate sectors
    parse(p,varargin{:});

    % Plot?
    plotLims=p.Results.Plot;
    if ~isempty(plotLims)
        figure
    end

    % Intialize output
    zSC=(NaN+NaN*1i)*zeros(size(uK));
    aSC=(NaN+NaN*1i)*zeros(size(uK));
    sigA=(NaN+NaN*1i)*zeros(size(uK));
    sigZ=(NaN+NaN*1i)*zeros(size(uK));
    rms=(NaN)*zeros(size(uK));

    % Loop over k-values
    for ik=1:length(uK)
        % Select energy vector and quasiparticle green's function
        zV=Z(ik,:);
        gV=QG(ik,:);

        % Remove NaN-values
        gV=gV(~isnan(zV));
        zV=zV(~isnan(zV));
        
        % Separate sectors
        if p.Results.SepReg
            if uK(ik) < kTh(1) % Sector B before energy minimum
                pts=(zV > th13(ik));
                cmplx=true;
            elseif uK(ik) < kTh(2) % Sector A before energy minimum
                pts=(th13(ik) > zV & zV > th12(ik));
                cmplx=true;
            elseif uK(ik) < kTh(3)
                % Real input values
                aIn=real(aIn);
                zIn=real(zIn);

                pts=(th12(ik) > zV);
                cmplx=false;
            elseif uK(ik) < kTh(4) % Sector A after energy minimum
                thTMP=[th13(ik),th12(ik)];
                pts=(max(thTMP) > zV & zV > min(thTMP));
                cmplx=true;
            else  % Sector B after energy minimum
                pts=(zV > th13(ik));
                cmplx=true;
            end
        else
            if kTh(2) < uK(ik) && uK(ik) < kTh(3)
                % Real input values
                aIn=real(aIn);
                zIn=real(zIn);
                pts=(th12(ik) > zV);
                cmplx=false;
            else
                pts=true(size(zV));
                cmplx=true;
            end
        end

        % If no valid data points were found, go to next k-value
        if ~any(pts)
            warning(['No valid points found for k=',num2str(uK(ik)),', returning NaN value.']);
            aSC(ik)=NaN+NaN*1i;
            zSC(ik)=NaN+NaN*1i;
            rms(ik)=NaN;
            sigA(ik)=NaN;
            sigZ(ik)=NaN;
            continue
        end

        % Fit
        if cmplx 
            [zFit,aFit,resnorm,~,residual,jac] ...
                =fitInvDet(zV(pts),gV(pts),zIn,aIn,p.Results.Filt);
        else
            zSCtmp=interp1(SC.kV,SC.zV,uK(ik),'makima');
            [zFit,aFit,resnorm,~,residual,jac] ...
                =fitInvDetRe(zV(pts),real(gV(pts)).^(-1),zIn,aIn,zSCtmp,p.Results.FiltRe);
        end

        % Use fit results for next starting values
        zIn=zFit;
        aIn=aFit;

        % Save final values
        aSC(ik,1)=aFit;
        zSC(ik,1)=zFit;
        rms(ik,1)=sqrt(resnorm/length(zV));

        % Calculate variance on A and Z
        xCovariance = inv(jac.'*jac)*var(residual);
        sigA(ik,1)=sqrt(xCovariance(1,1));
        sigZ(ik,1)=sqrt(xCovariance(2,2));

        % Plot results if asked for
        if ~isempty(plotLims) && plotLims(1) < uK(ik) && uK(ik) < plotLims(2)
                
            % real and imaginary part above th12
            if uK(ik) < kTh(2) || kTh(3) < uK(ik)

                % Calculate fitting region and make fit vector for plot
                [yMax,~]=max(imag(gV(pts)));
                iData=find(imag(gV(pts))>=yMax*p.Results.Filt);
                zVPts=zV(pts);
                zVFit=linspace(zVPts(iData(1)),zVPts(iData(end)),101);

                fitFun=aFit./(zFit-zVFit);
                fitPtRe= imag(aFit)/imag(zFit);
                fitPtIm=-real(aFit)/imag(zFit);

                % Data points
                zDatEndPts=[zVPts(iData(1)),zVPts(iData(end))];

                subplot(1,2,1)
                    plot(zVFit,real(fitFun) ...
                        ,zV,real(gV),'.' ...
                        ,real(zFit),fitPtRe,'o' ...
                    ,'LineWidth',1.5);
                    hold on
                    plotKline([th12(ik),th13(ik)], ...
                        [min(real(gV)),max(real(gV))]);
                    % plot(zDatEndPts,[0,0],'LineWidth',5);
                    hold off
                    title(['th12=',num2str(th12(ik))]);
                subplot(1,2,2)
                    plot(zVFit,imag(fitFun) ...
                        ,zV,imag(gV),'.' ...
                        ,real(zFit),fitPtIm,'o' ...
                    ,'LineWidth',1.5);
                    hold on
                    plotKline([th12(ik),th13(ik)], ...
                        [0,max(imag(gV))]);
                    % plot(zDatEndPts,[0,0],'LineWidth',5);
                    hold off
                title(['k=',num2str(uK(ik))]);
                pause
            else %Only real part

                % Calculate fitting region and make fit vector for plot
                i1=find(zSCtmp-p.Results.FiltRe<zV(pts),1,'first');
                i2=find(zV(pts)<zSCtmp+p.Results.FiltRe,1,'last');
                iData=[i1,i2];
                zVPts=zV(pts);
                zVFit=linspace(zVPts(iData(1)),zVPts(iData(end)),101);

                plot(zVFit,real((zFit-zVFit)./aFit) ...
                    ,zV,real(gV).^(-1),'.','LineWidth',1.5);
                title(['k=',num2str(uK(ik))])
                hold on
                plotKline(th12(ik),[-.1,.1]);
                hold off
                pause
            end
        end
    
    end

end

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

function [zFit,aFit,resnorm,output,residual,jac]=fitInvDet(zV,gV,zIn,aIn,fltr)
% Fit the inverse of the determinant
%   gV = aFit/(zV-zFit)
% with aFit and zFit complex numbers.
% zIn and aIn are the complex initial values

    % Fit options
    opts = optimoptions(@lsqcurvefit,'Display','off');

    % Function to fit
    fitFun=@(v,epsilon) v(1)./(v(2)-epsilon);

    % Initial guess
    v0=[aIn;zIn];

    % Only look at points above yMax*fltr
    yMax=max(imag(gV));
    iData=find(imag(gV)>=yMax*fltr);
    xData=zV(iData);
    yData=gV(iData);

    % Fit function
    [vestimated,resnorm,residual,~,output,~,jac] = lsqcurvefit(fitFun,v0,xData,yData,[],[],opts);

    % Output
    aFit=vestimated(1);
    zFit=vestimated(2);

end


function [zFit,aFit,resnorm,output,residual,jac]=fitInvDetRe(zV,gV,zIn,aIn,zSC,fltr)
% Fit the inverse of the determinant
%   gV = aFit/(zV-zFit)
% with aFit and zFit real numbers.
% zIn and aIn are the real initial values

    % Fit options
    opts = optimoptions(@lsqcurvefit,'Display','off');

    % Function to fit
    fitFun=@(v,epsilon) (v(2)-epsilon)./v(1);

    % Initial guess
    v0=real([aIn;zIn]);

    % Only look at points close to zero
    iData=(zSC-fltr<zV & zV<zSC+fltr);
    xData=zV(iData);
    yData=gV(iData);

    % Fit function
    [vestimated,resnorm,residual,~,output,~,jac] = lsqcurvefit(fitFun,v0,xData,yData,[],[],opts);

    % Output
    aFit=vestimated(1);
    zFit=vestimated(2);

end