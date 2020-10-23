%% MANAGE workspaces
close all; % close all figures
clear all; % clear the workspace
clc;

% Please, change the paths !
eval('pathManager');

multicpu = true;

%% READ parfile
parFileName = 'parFileLOOPS';
eval(parFileName);

%% INSTANTIATING THE OOMAO CLASSES

% source
ngs = source('wavelength',photometry.V0);

% telescope
tel = telescope(D,'obstructionRatio',cobs,'resolution',nPxPup);
pupZeroPad = padarray(tel.pupil,[(Samp+1)*nPxPup/2,(Samp+1)*nPxPup/2]);

atm = atmosphere(photometry.V0,r0_mean,L0,'altitude',0,'fractionnalR0',1,...
    'windSpeed',v_mean,'windDirection',0);

% wavefront-sensor
pyr = pyramid(nLenslet,nPxPup,'modulation',pyrMod,'binning',pyrBinning,'c',Samp);
ngs = ngs.*tel*pyr;
pyr.INIT;

% spatial frequency model
fao = spatialFrequencyAdaptiveOptics(tel,atm,nLenslet+1,noiseVar,...
    loopGain, sampTime, latency, resAO,'pyramid',0,'modulation',pyrMod,...
    'nTimes',nTimes);
            
% AO correction area deifnition
[fxExt,fyExt] = freqspace(size(fao.fx,1)*fao.nTimes,'meshgrid');
fxExt = fxExt*fao.fc*fao.nTimes;
fyExt = fyExt*fao.fc*fao.nTimes;
index = abs(fxExt)<fao.fc & abs(fyExt)<fao.fc;

% Fitting, Noise and Kolmogorov PSDs
psdFit = fao.fittingPSD(fao.fx,fao.fy)/fao.atm.r0^(-5/3);
psdNoise = fao.noisePSD(fao.fx,fao.fy)/fao.noiseVariance;

%mfew parameters to create the phase screen from the PSD
N       = 4*nPxPup;
L       = (N-1)*tel.D/(nPxPup-1);
[fx,fy] = freqspace(N,'meshgrid');
[~,fr]  = cart2pol(fx,fy);
fr      = fftshift(fr.*(N-1)/L./2);
[idx]           = find(fr==0);
fourierSampling = 1./L;

fc  = fao.fc;
mskIn = ~(abs(fao.fx)>=fc | abs(fao.fy)>=fc);
Ts  = fao.samplingTime;
fxx = fao.fx;
fyy = fao.fy;
F = (fao.Rx.*fao.SxAv + fao.Ry.*fao.SyAv);
psdKolmo = fao.pistonFilter(hypot(fxx,fyy)).*phaseStats.spectrum(hypot(fxx,fyy),fao.atm)/fao.atm.r0^(-5/3);

close all;
%% CREATING LINEAR HYPERBOLIC SAMPLING FOR INPUTS PARAMETERS - MAXIMIN - TRUNCATED 

r0_list = lhsdesign_modified(nr0,r0_min,r0_mean + r0_std);
v_list  = lhsdesign_modified(nv,v_min,v_mean + v_std);
n_list  = lhsdesign_modified(nn,n_min,n_mean + n_std);
nAll    = nr0*nv*nn

%% LOOP - SINGLE CPU

if ~multicpu
    kIter = 0;
    elapsed = 0;
    
    for iv = 1:nv
        % update the windspeed value
        atm.layer.windSpeed = v_list(iv);
        % update the AO controller model
        fao.controller('int');
        for ir0 = 1:1
            %update the r0 value
            atm.r0 = r0_list(ir0);
            psdServo = fao.anisoServoLagPSD(fao.fx,fao.fy);
            for in = 1:1
                kIter = kIter + 1;
                t_i = tic();
                
                % Get the residual phase PSD
                psdAO = psdFit*atm.r0^(-5/3);
                psdAO(index) = psdAO(index) + reshape(psdNoise*n_list(in) + psdServo,[],1);
                
                % Get the phase map
                map = real(ifft2(idx.*sqrt(fftshift(psdAO)).*fft2(randn(atm.rngStream,N))./N).*fourierSampling).*N.^2;
                phaseMap = map(1:nPxPup,1:nPxPup);
                
                % propagate through the pyramid
                ngs = ngs.*tel;
                ngs.phase = phaseMap;
                ngs = ngs*pyr;
                
                % record outputs
                fitswrite(single(phaseMap),[path_save,'ground_truth_phase_r0_',num2str(r0_list(ir0)),'_v_',num2str(v_list(iv)),'_noise_',num2str(n_list(in)),'.fits']);
                fitswrite(single(pyr.camera.frame),[path_save,'measurements_intensity_r0_',num2str(r0_list(ir0)),'_v_',num2str(v_list(iv)),'_noise_',num2str(n_list(in)),'.fits']);
                elapsed = elapsed + toc(t_i);
                fprintf(['Time remaining :',num2str(elapsed/kIter*(nAll-kIter)/3600),' h\n']);
                
            end
        end
    end

    %
    close all;
    
    figure;
    imagesc(ngs.phase);
    pbaspect([1,1,1]);
    colorbar('TickLabelInterpreter','latex')
    xlabel('Pixel in the pupil','interpreter','latex','fontsize',16)
    ylabel('Pixel in the pupil','interpreter','latex','fontsize',16)
    set(gca,'FontSize',16,'FontName','cmr12','TickLabelInterpreter','latex');
    
    figure;
    imagesc(pyr.camera.frame);
    pbaspect([1,1,1]);
    colorbar('TickLabelInterpreter','latex')
    xlabel('Pixel in the WFS detector plane','interpreter','latex','fontsize',16)
    ylabel('Pixel in the WFS detector plane','interpreter','latex','fontsize',16)
    set(gca,'FontSize',16,'FontName','cmr12','TickLabelInterpreter','latex');


else
    %% LOOP - MULTI CPU
    
    % Check the presence of the COntrol System tool box
    toolboxesList = ver;
    flagControlToolBox = any(strcmp(cellstr(char(toolboxesList.Name)), 'Control System Toolbox'));
    
    % Extending list of parameters to span a single vector
    r0_list_ext = repelem(r0_list,nv*nn);
    v_list_ext = repelem(repmat(v_list,nr0,1),nn);
    n_list_ext = repmat(n_list,nr0*nv,1);
    
    % Initialization of system parameters
    nPts        = size(fao.fx,1);
    nTh_        = fao.nTh;
    h1buf       = zeros(nPts,nPts,nTh_);
    h2buf       = zeros(nPts,nPts,nTh_);
    thetaWind   = linspace(0, 2*pi-2*pi/nTh_,nTh_);
    costh       = cos(thetaWind);    
    RTF         = fao.atf;
    rngStream   = atm.rngStream;
    pupil       = tel.pupil;
    
    % Give workers access to OOMAO functions
    addAttachedFiles(gcp,{'telescope.m','telescopeAbstract.m','pyramid.m','source.m'})
    
    t1 = tic();
    
    parfor kIter = 1:10
        
        t_i = tic();
        
        % Get the fitting PSD
        r053  = r0_list_ext(kIter)^(-5/3);
        psdAO = psdFit*r053;
        
        % Get the noise PSD
        psdN = psdNoise*n_list_ext(kIter);
        
        % Get the servo-lag PSDs
        h1buf  = zeros(nPts,nPts,nTh_);
        h2buf = zeros(nPts,nPts,nTh_);
        for iTheta = 1:nTh_
            fi = -v_list_ext(kIter)*fxx*costh(iTheta);
            id = abs(fi) <1e-7;
            fi(id) = 1e-8.*sign(fi(id));
            if flagControlToolBox
                [MAG, PH] = bode(RTF, 2*pi*fi(:));
                MAG = reshape(MAG,[nPts,nPts]);
                MAG(fi == 0) = 1;
                PH = reshape(PH,[nPts,nPts]);
                h2buf(:,:,iTheta) = abs(MAG).^2;
                h1buf(:,:,iTheta) = MAG.*exp(1i*PH/180*pi);
            else
                tmp = RTF(exp(-2i*pi*fi*Ts));
                h1buf(:,:,iTheta) = tmp;
                h2buf(:,:,iTheta) = abs(tmp).^2;
            end
        end
        h1 = sum(h1buf,3)/nTh_;
        h2 = sum(h2buf,3)/nTh_;
        
        psdServo = (1 + abs(F).^2.*h2 - F.*h1 - conj(F.*h1)).*psdKolmo;
        
        % Get the residual phase PSD
        psdAO(index) = psdAO(index) + reshape(psdN(mskIn) + psdServo(mskIn),[],1);
        
        % Get the phase map
        map = real(ifft2(idx.*sqrt(fftshift(psdAO)).*fft2(randn(rngStream,N))./N).*fourierSampling).*N.^2;
        phaseMap = pupil.*map(1:nPxPup,1:nPxPup);
        
        % propagate through the pyramid
        n2 = times(ngs,tel);
        n2.phase = phaseMap;
        n2 = mtimes(n2,pyr);
        
        % record outputs
        fitswrite(single(phaseMap),[path_save,'ground_truth_phase_r0_',num2str(r0_list_ext(kIter)),'_v_',num2str(v_list_ext(kIter)),'_noise_',num2str(n_list_ext(kIter)),'.fits']);
        fitswrite(single(pyr.camera.frame),[path_save,'measurements_intensity_r0_',num2str(r0_list_ext(kIter)),'_v_',num2str(v_list_ext(kIter)),'_noise_',num2str(n_list_ext(kIter)),'.fits']);
        elapsed = toc(t_i);
        fprintf(['Time remaining :',num2str(elapsed*nAll/kIter/3600),' h\n']);
    end
    tt = toc(t1);
end


