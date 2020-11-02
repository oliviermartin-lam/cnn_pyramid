%% MANAGE workspaces
close all; % close all figures
clear all; % clear the workspace
clc;

% Please, change the paths !
eval('pathManager');

%% READ parfile
parFileName = 'parFileLOOPS';
eval(parFileName);

%% INSTANTIATING THE OOMAO CLASSES

% source
ngs = source('wavelength',photoNGS);

% telescope
tel = telescope(D,'obstructionRatio',cobs,'resolution',nPxPup);
pupZeroPad = padarray(tel.pupil,[(Samp+1)*nPxPup/2,(Samp+1)*nPxPup/2]);

atm = atmosphere(photoNGS,(r0_max-r0_min)/2,L0,'altitude',0,'fractionnalR0',1,...
    'windSpeed',(v_max - v_min)/2,'windDirection',0);

% wavefront-sensor
pyr = pyramid(nLenslet,nPxPup,'modulation',pyrMod,'binning',pyrBinning,'c',Samp);
% Flat for FULL FRAME
ngs = ngs.*tel*pyr;
pyr.INIT
ngs = ngs.*tel*pyr;
I_0 = pyr.camera.frame./sum(pyr.camera.frame(:));

% spatial frequency model
fao = spatialFrequencyAdaptiveOptics(tel,atm,nLenslet+1,(n_max-n_min)/2,...
    loopGain, sampTime, latency, resAO,'pyramid',0,'modulation',pyrMod,...
    'nTimes',nTimes);
close all;
            
% AO correction area definition
[fxExt,fyExt] = freqspace(size(fao.fx,1)*fao.nTimes,'meshgrid');
fxExt = fxExt*fao.fc*fao.nTimes;
fyExt = fyExt*fao.fc*fao.nTimes;
index = abs(fxExt)<fao.fc & abs(fyExt)<fao.fc;

% Fitting, Noise and Kolmogorov PSDs
psdFit = fao.fittingPSD(fao.fx,fao.fy)/fao.atm.r0^(-5/3);
psdNoise = fao.noisePSD(fao.fx,fao.fy)/fao.noiseVariance;

%mfew parameters to create the phase screen from the PSD
N       = 2*Samp*nPxPup;
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


%% ZERNIKE RECONSTRUCTION MATRIX : FULL-FRAME TO ZERNIKE MODE
% define Zernike Modes
zern     = zernike(2:nZern+1,tel.D,'resolution',tel.resolution);
iMat     = interaction_matrix(ngs,tel,pyr,zern.modes);
pyr2zern = pinv(iMat);
  
    
% slopes to zernike
%dm      = deformableMirror(nZern,'modes',zern,'resolution',tel.resolution,'validActuator',true(1,nZern));
%ngs     = ngs.*tel;
%calibDm = calibration(dm,pyr,ngs,ngs.wavelength/1000);
%pyr2zern= calibDm.M;

%% CREATING LINEAR HYPERBOLIC SAMPLING FOR INPUTS PARAMETERS - MAXIMIN - TRUNCATED 
nAll     = 1250000;
r0_list  = r0_min + (r0_max - r0_min)*rand(1,nAll);
v_list   = v_min + (v_max - v_min)*rand(1,nAll);
n_list   = n_min + (n_max - n_min)*rand(1,nAll);

%% CHECK THE RESIDUAL PHASE VARIANCE AT MEAN CONDITIONS

% mean conditions
fao.atm.r0 = mean(r0_list);
psdAO_mean = zeros(size(fxExt));
psdAO_mean(index) = fao.anisoServoLagPSD()  + psdNoise*mean(n_list);            
psdAO_mean = psdAO_mean + psdFit*mean(r0_list)^(-5/3);
            
nIter = 1e3;
wfe   = zeros(1,nIter);
for k=1:nIter
    % mean
    phaseMap = real(ifft2(idx.*sqrt(fftshift(psdAO_mean)).*fft2(randn(atm.rngStream,N))./N).*fourierSampling).*N.^2;
    phaseMap = tel.pupil.*phaseMap(1:nPxPup,1:nPxPup);
    wfe(k)   = std(phaseMap(tel.pupilLogical(:)))*atm.wavelength*1e9/2/pi;
end

figure;
histogram(wfe,'Normalization','Probability','BinEdges', min(wfe):10:max(wfe));hold on;
ylabel('Probability','interpreter','latex','fontsize',20);
xlabel('Wavefront error (nm)','interpreter','latex','fontsize',20);
set(gca,'FontSize',20,'FontName','cmr12','TickLabelInterpreter','latex');
pbaspect([1.6,1,1]);


%% CHECK THE WFS MEASUREMENTS
phaseMap = real(ifft2(idx.*sqrt(fftshift(psdAO_mean)).*fft2(randn(atm.rngStream,N))./N).*fourierSampling).*N.^2;
phaseMap = tel.pupil.*phaseMap(1:nPxPup,1:nPxPup);
ngs       = ngs.*tel;
ngs.phase = phaseMap;
ngs       = ngs*pyr;

figure;
subplot(1,2,1)
imagesc(ngs.phase);
title('Input phase in rad','interpreter','latex');
pbaspect([1,1,1]);
colorbar('TickLabelInterpreter','latex')
xlabel('Pixel in the pupil','interpreter','latex','fontsize',12)
ylabel('Pixel in the pupil','interpreter','latex','fontsize',12)
set(gca,'FontSize',12,'FontName','cmr12','TickLabelInterpreter','latex');

subplot(1,2,2)
imagesc(pyr.camera.frame/sum(pyr.camera.frame(:)) - I_0);
title('Pyramid measurements ','interpreter','latex');
pbaspect([1,1,1]);
colorbar('TickLabelInterpreter','latex')
xlabel('Pixel in the WFS detector plane','interpreter','latex','fontsize',12)
ylabel('Pixel in the WFS detector plane','interpreter','latex','fontsize',12)
set(gca,'FontSize',12,'FontName','cmr12','TickLabelInterpreter','latex');

% subplot(1,3,3)
% imagesc(int8(pyr.camera.frame/sum(pyr.camera.frame(:)) - I_0));
% title('8-bits pyramid measurements','interpreter','latex');
% pbaspect([1,1,1]);
% colorbar('TickLabelInterpreter','latex')
% xlabel('Pixel in the WFS detector plane','interpreter','latex','fontsize',12)
% ylabel('Pixel in the WFS detector plane','interpreter','latex','fontsize',12)
% set(gca,'FontSize',12,'FontName','cmr12','TickLabelInterpreter','latex');

%% CHECK ZERNIKE RECONSTRUCTION
close all;clc;

% if you change the phase amplitude, you'll see that the pyramid-based
% Zernike reconstruction degrades -> optical gain issue
wvl_factor  = ngs.wavelength*1e9/2/pi; % from rad 2 nm
amp_z       = 1/wvl_factor;
z_true      = randn(nZern,1)*amp_z;
ngs         = ngs.*tel;
ngs.phase   = reshape(zern.modes*z_true,tel.resolution,[]);
ngs         = ngs*pyr;

% phase map to zernike
ph2zern     = pinv(zern.modes);
c1          = wvl_factor*ph2zern*ngs.phase(:);
zMap1       = reshape(sum(zern.modes*c1,2),tel.resolution,[]);
zMap1       = zMap1 - mean(zMap1(tel.pupilLogical(:)));

% pyramid full-frame to zernike
pyr_frame   = pyr.camera.frame./sum(pyr.camera.frame(:))-I_0;
c2          = wvl_factor*pyr2zern*pyr_frame(:);
zMap2       = reshape(sum(zern.modes*c2,2),tel.resolution,[]);
zMap2       = zMap2 - mean(zMap2(tel.pupilLogical(:)));


% check numbers
sqrt(ngs.var)*wvl_factor
sqrt(sum(c1.^2))
sqrt(sum(c2.^2))

figure;
plot(2:nZern+1,z_true*wvl_factor,'k-','linewidth',1.5);
hold on;
plot(2:nZern+1,c1,'b--','linewidth',1.5);
plot(2:nZern+1,c2,'r--','linewidth',1.5);
xlabel('Noll''s j-index','interpreter','latex','fontsize',20);
ylabel('Zernike coefficients (nm)','interpreter','latex','fontsize',20);
set(gca,'FontSize',20,'FontName','cmr12','TickLabelInterpreter','latex');
legend({'Ground truth','Phase map-based reconstruction','Pyramid-slopes-based reconstruction'},'interpreter','latex','fontsize',16);


figure;
subplot(2,3,1)
imagesc(1e9*ngs.meanRmOpd);
title('Input OPD in nm','interpreter','latex');
pbaspect([1,1,1]);
colorbar('TickLabelInterpreter','latex')
xlabel('Pixel in the pupil','interpreter','latex','fontsize',12)
ylabel('Pixel in the pupil','interpreter','latex','fontsize',12)
set(gca,'FontSize',12,'FontName','cmr12','TickLabelInterpreter','latex');

subplot(2,3,2)
imagesc(zMap1);
title('Zernike reconstruction from the phase ','interpreter','latex');
pbaspect([1,1,1]);
colorbar('TickLabelInterpreter','latex')
xlabel('Pixel in the pupil','interpreter','latex','fontsize',12)
ylabel('Pixel in the pupil','interpreter','latex','fontsize',12)
set(gca,'FontSize',12,'FontName','cmr12','TickLabelInterpreter','latex');

subplot(2,3,3)
imagesc(zMap2);
title('Zernike reconstruction from the pyramid ','interpreter','latex');
pbaspect([1,1,1]);
colorbar('TickLabelInterpreter','latex')
xlabel('Pixel in the pupil','interpreter','latex','fontsize',12)
ylabel('Pixel in the pupil','interpreter','latex','fontsize',12)
set(gca,'FontSize',12,'FontName','cmr12','TickLabelInterpreter','latex');

subplot(2,3,5)
imagesc(zMap1 - 1e9*ngs.meanRmOpd);
title('Residual','interpreter','latex');
pbaspect([1,1,1]);
colorbar('TickLabelInterpreter','latex')
xlabel('Pixel in the pupil','interpreter','latex','fontsize',12)
ylabel('Pixel in the pupil','interpreter','latex','fontsize',12)
set(gca,'FontSize',12,'FontName','cmr12','TickLabelInterpreter','latex');

subplot(2,3,6)
imagesc(zMap2 - 1e9*ngs.meanRmOpd);
title('Residual','interpreter','latex');
pbaspect([1,1,1]);
colorbar('TickLabelInterpreter','latex')
xlabel('Pixel in the pupil','interpreter','latex','fontsize',12)
ylabel('Pixel in the pupil','interpreter','latex','fontsize',12)
set(gca,'FontSize',12,'FontName','cmr12','TickLabelInterpreter','latex');

%% LOOP - MULTI CPU

% Check the presence of the COntrol System tool box
toolboxesList = ver;
flagControlToolBox = any(strcmp(cellstr(char(toolboxesList.Name)), 'Control System Toolbox'));

% Initialization of system parameters
nPts        = size(fao.fx,1);
RTF         = fao.atf;
pupil       = tel.pupil;
rngStream   = atm.rngStream;

% indexes for cropping the pyramid camera frame
idx1 = ((pyr.c-1)/2) * pyr.nLenslet + 1 : ((pyr.c-1)/2 + 1) * pyr.nLenslet;
idx2 = ((pyr.c-1)/2 + pyr.c) * pyr.nLenslet + 1 : ((pyr.c-1)/2 + pyr.c + 1) * pyr.nLenslet;

% Give workers access to OOMAO functions
addAttachedFiles(gcp,{'telescope.m','telescopeAbstract.m','pyramid.m','source.m'})
t1 = tic();


parfor kIter = 1:nAll
    
    t_i = tic();
    
    % Get the fitting PSD
    psdAO = psdFit*r0_list(kIter)^(-5/3);
    
    % Get the noise PSD
    psdN = psdNoise*n_list(kIter);
    
    % Get the servo-lag PSDs
    fi      = -v_list(kIter)*fxx;
    id      = abs(fi) <1e-7;
    fi(id)  = 1e-8.*sign(fi(id));
    
    if flagControlToolBox
        [MAG, PH] = bode(RTF, 2*pi*fi(:));
        MAG = reshape(MAG,[nPts,nPts]);
        MAG(fi == 0) = 1;
        PH = reshape(PH,[nPts,nPts]);
        h2 = abs(MAG).^2;
        h1 = MAG.*exp(1i*PH/180*pi);
    else
        h1 = RTF(exp(-2i*pi*fi*Ts));
        h2 = abs(h1).^2;
    end
    
    psdServo = (1 + abs(F).^2.*h2 - F.*h1 - conj(F.*h1)).*psdKolmo;
    
    % Get the residual phase PSD
    psdAO(index) = psdAO(index) + reshape(psdN(mskIn) + psdServo(mskIn),[],1);
    
    % Get the phase map
    map = real(ifft2(idx.*sqrt(fftshift(psdAO)).*fft2(randn(rngStream,N))./N).*fourierSampling).*N.^2;
    phaseMap = pupil.*map(1:nPxPup,1:nPxPup);
    
    % propagate through the pyramid
    n2          = times(ngs,tel);
    n2.phase    = phaseMap;
    n2          = mtimes(n2,pyr);
    pyr_frame   = pyr.camera.frame./sum(pyr.camera.frame(:))-I_0;
    
    % Ground-truth Zernike reconstruction in nm
    zCoefs = wvl_factor*ph2zern*n2.phase(:);

    % Pyramid-based Zernike reconstruction in nm
    zCoefs_pyr =  wvl_factor*pyr2zern*pyr_frame(:);
    
    % wavefront-error
    wfe = sqrt(n2.var)*wvl_factor;

    % crop image
    pyr_frame = ([pyr_frame(idx1,idx1),pyr_frame(idx2,idx1);pyr_frame(idx1,idx2),pyr_frame(idx2,idx2)]);
    
    % record outputs
    fitswrite(single(phaseMap),[path_save,'ground_truth_phase_r0_',num2str(r0_list(kIter)),'_v_',num2str(v_list(kIter)),'_noise_',num2str(n_list(kIter)),'_wfe_',num2str(wfe),'nm.fits']);
    fitswrite(single(zCoefs),[path_save,'ground_truth_zernike_r0_',num2str(r0_list(kIter)),'_v_',num2str(v_list(kIter)),'_noise_',num2str(n_list(kIter)),'_wfe_',num2str(wfe),'nm.fits']);
    fitswrite(single(zCoefs_pyr),[path_save,'pyramid_based_zernike_r0_',num2str(r0_list(kIter)),'_v_',num2str(v_list(kIter)),'_noise_',num2str(n_list(kIter)),'_wfe_',num2str(wfe),'nm.fits']);
    fitswrite(single(pyr_frame),[path_save,'measurements_intensity_r0_',num2str(r0_list(kIter)),'_v_',num2str(v_list(kIter)),'_noise_',num2str(n_list(kIter)),'_wfe_',num2str(wfe),'nm.fits']);
    elapsed = toc(t_i);
    fprintf(['Time remaining :',num2str(elapsed*nAll/kIter/3600),' h\n']);
end
tt = toc(t1);


