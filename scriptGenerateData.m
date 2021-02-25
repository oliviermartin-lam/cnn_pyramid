%% MANAGE workspaces
close all; % close all figures
clear all; % clear the workspace
clc;

% Please, change the paths !
eval('pathManager');

%% READ parfile
parFileName = 'parFileLOOPS';
eval(parFileName);
check = false;

%% INSTANTIATING THE OOMAO CLASSES

% source
ngs = source('wavelength',photoNGS,'magnitude',magNGS);

% telescope
tel = telescope(D,'obstructionRatio',cobs,'resolution',nPxPup);
pupZeroPad = padarray(tel.pupil,[(Samp+1)*nPxPup/2,(Samp+1)*nPxPup/2]);

atm = atmosphere(photoNGS,r0,L0,'altitude',0,'fractionnalR0',1,...
    'windSpeed',wSpeed,'windDirection',wDirection*pi/180);

% wavefront-sensor
pyr = pyramid(nLenslet,nPxPup,'modulation',pyrMod,'binning',pyrBinning,'c',Samp);
% Flat for FULL FRAME
ngs = ngs.*tel*pyr;
pyr.INIT
ngs = ngs.*tel*pyr;
I_0 = pyr.camera.frame./sum(pyr.camera.frame(:));

wvl = ngs.wavelength;

%% Phase screen generator

if strcmpi(genType,'FOURIER')
    fao = spatialFrequencyAdaptiveOptics(tel,atm,nLenslet+1,noiseVar,...
        loopGain, sampTime, latency, resAO,'pyramid',0,'modulation',pyrMod,...
        'nTimes',nTimes);
    close all;
    
    [fxExt,fyExt] = freqspace(size(fao.fx,1)*fao.nTimes,'meshgrid');
    fxExt = fxExt*fao.fc*fao.nTimes;
    fyExt = fyExt*fao.fc*fao.nTimes;
    index = abs(fxExt)<fao.fc & abs(fyExt)<fao.fc;
    psdKolmo = fao.pistonFilter(hypot(fxExt,fyExt)).*phaseStats.spectrum(hypot(fxExt,fyExt),fao.atm);
    
    if loopGain == 0
        psd = psdKolmo;
        
    else
        psd        = zeros(size(fxExt));
        psd(index) = fao.noisePSD(fao.fx,fao.fy) + fao.anisoServoLagPSD(fao.fx,fao.fy);
        psd        = psd + fao.fittingPSD(fao.fx,fao.fy);
        
    end
    % Zernike modes to be reconstructed
    jIndex = 2:nZern+1;
   
    % few parameters to create the phase screen from the PSD
    N       = 2*Samp*nPxPup;
    L       = (N-1)*tel.D/(nPxPup-1);
    [fx,fy] = freqspace(N,'meshgrid');
    [~,fr]  = cart2pol(fx,fy);
    fr      = fftshift(fr.*(N-1)/L./2);
    [idx]           = find(fr==0);
    fourierSampling = 1./L;
    
    % Check the presence of the Control System tool box
    toolboxesList = ver;
    flagControlToolBox = any(strcmp(cellstr(char(toolboxesList.Name)), 'Control System Toolbox'));
    
    % Initialization of system parameters
    nPts        = size(fao.fx,1);
    RTF         = fao.atf;
    pupil       = tel.pupil;
    rngStream   = atm.rngStream;
    
else
    % ZERNIKE CLASS
    zernGen = zernike(jIndex,tel.D,'resolution',tel.resolution);
    nZern   = numel(jIndex);
    % ZERNIKE POLYNOMES  - ARRAY OF SIZE nPix^2 x nModes
    zModes  = zernGen.modes;
    % DISTRIBUTION OF AMPLITUDE
    if numel(zStdMax) == 1
        zStdMax = zStdMax * ones(nZern,1);
    else
        zStdMax = reshape(zStdMax,nZern,1);
    end
    if numel(zMean) == 1
        zMean = zMean * ones(nZern,1);
    else
        zMean = reshape(zMean,nZern,1);
    end
    
    if strcmpi(zDistrib,'NORMAL')
        % Normal distribution of each of the nZ_Gen modes
        zAmplitude =  bsxfun(@plus,zMean,bsxfun(@times,zStdMax,randn(nZern,nAll)));
    elseif strcmpi(zDistrib,'UNIFORM')
        zAmplitude = zMean - 5*zStdMax + 2*(5*zStdMax - zMean)*rand(nZern,nAll);
    end
end

%% ZERNIKE RECONSTRUCTION MATRIX : FULL-FRAME TO ZERNIKE MODE

% define Zernike Modes
zernRec  = zernike(jIndex,tel.D,'resolution',tel.resolution);
iMat     = interaction_matrix(ngs,tel,pyr,zernRec.modes);
pyr2zern = pinv(iMat);
wvl_factor  = wvl*1e9/2/pi; % from rad 2 nm
ph2zern     = pinv(zernRec.modes);

% if you change the phase amplitude, you'll see that the pyramid-based
% Zernike reconstruction degrades -> optical gain issue
if check
    close all;
    amp_z       = 1000/wvl_factor; %50 nm in amplitude
    z_true      = randn(nZern,1)*amp_z;
    ngs         = ngs.*tel;
    ngs.phase   = reshape(zernRec.modes*z_true,tel.resolution,[]);
    ngs         = ngs*pyr;
    
    % phase map to zernike
    c1          = wvl_factor*ph2zern*ngs.phase(:);
    zMap1       = reshape(sum(zernRec.modes*c1,2),tel.resolution,[]);
    zMap1       = zMap1 - mean(zMap1(tel.pupilLogical(:)));
    
    % pyramid full-frame to zernike
    pyr_frame   = pyr.camera.frame./sum(pyr.camera.frame(:))-I_0;
    c2          = wvl_factor*pyr2zern*pyr_frame(:);
    zMap2       = reshape(sum(zernRec.modes*c2,2),tel.resolution,[]);
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
end

%% LOOP - MULTI CPU
% indexes for cropping the pyramid camera frame
idx1 = ((pyr.c-1)/2) * pyr.nLenslet + 1 : ((pyr.c-1)/2 + 1) * pyr.nLenslet;
idx2 = ((pyr.c-1)/2 + pyr.c) * pyr.nLenslet + 1 : ((pyr.c-1)/2 + pyr.c + 1) * pyr.nLenslet;

% Give workers access to OOMAO functions
addAttachedFiles(gcp,{'telescope.m','telescopeAbstract.m','pyramid.m','source.m'})
t1 = tic();

% Manage noise
if ron
    pyr.camera.readOutNoise = ron;
    pyr.camera.photonNoise = true;
    pyr.camera.quantumEfficiency = QE;
end


if strcmpi(genType,'ZERNIKE')
   
    parfor kIter = 1:nAll
        
        t_i = tic();
        
        % ------------------- ZERNIKE GENERATOR
        zCoefs   = zAmplitude(:,kIter);
        phaseMap = 2*pi*1e-9/wvl*reshape(zModes*zCoefs,nPxPup,nPxPup);
        %id
        idk = ['jIndex_',sprintf('%d_', jIndex),'zCoefs_',sprintf('%.1f_', zCoefs)];
        
        
        % propagate through the pyramid
        n2          = times(ngs,tel);
        n2.phase    = phaseMap;
        n2          = mtimes(n2,pyr);
        pyr_frame   = pyr.camera.frame./sum(pyr.camera.frame(:))-I_0;
        
        
        % Ground-truth - FOURIER GENERATOR CASE
        if strcmpi(genType,'FOURIER')
            zCoefs = wvl_factor*ph2zern*n2.phase(:);
        end
        
        % Pyramid-based Zernike reconstruction in nm
        zCoefs_pyr =  wvl_factor*pyr2zern*pyr_frame(:);
        
        % wavefront-error
        wfe = sqrt(n2.var)*wvl_factor;
        
        % crop image
        pyr_frame = ([pyr_frame(idx1,idx1),pyr_frame(idx2,idx1);pyr_frame(idx1,idx2),pyr_frame(idx2,idx2)]);
        
        % record outputs
        % the phase map recording is disbaled to reduce the data amount to be saved
        %fitswrite(single(phaseMap),[path_save,'ground_truth_phase_',idk,'_wfe_',sprintf('%.1f_', wfe),'nm.fits']);
        fitswrite(single(zCoefs),[path_save,'ground_truth_zernike_',idk,'_wfe_',sprintf('%.1f_', wfe),'nm.fits']);
        fitswrite(single(zCoefs_pyr),[path_save,'pyramid_based_zernike_',idk,'_wfe_',sprintf('%.1f_', wfe),'nm.fits']);
        fitswrite(single(pyr_frame),[path_save,'measurements_intensity_',idk,'_wfe_',sprintf('%.1f_', wfe),'nm.fits']);
        elapsed = toc(t_i);
        fprintf(['Time remaining :',num2str(elapsed*nAll/kIter/3600),' h\n']);
    end
else
    
    parfor kIter = 1:nAll
        
        t_i = tic();
        
        % ------------------- FOURIER GENERATOR
        % Get the phase map
        map = real(ifft2(idx.*sqrt(fftshift(psd)).*fft2(randn(rngStream,N))./N).*fourierSampling).*N.^2;
        phaseMap = pupil.*map(1:nPxPup,1:nPxPup);
        %id
        idk = ['r0_',num2str(r0),'_v_',num2str(wSpeed),'_noise_',num2str(noiseVar)];
        
        
        % propagate through the pyramid
        n2          = times(ngs,tel);
        n2.phase    = phaseMap;
        n2          = mtimes(n2,pyr);
        pyr_frame   = pyr.camera.frame./sum(pyr.camera.frame(:))-I_0;    
        
        % Ground-truth - FOURIER GENERATOR CASE
        zCoefs = wvl_factor*ph2zern*n2.phase(:);
           
        % Pyramid-based Zernike reconstruction in nm
        zCoefs_pyr =  wvl_factor*pyr2zern*pyr_frame(:);
        
        % wavefront-error
        wfe = sqrt(n2.var)*wvl_factor;
        
        % crop image
        pyr_frame = ([pyr_frame(idx1,idx1),pyr_frame(idx2,idx1);pyr_frame(idx1,idx2),pyr_frame(idx2,idx2)]);
        
        % record outputs
        % the phase map recording is disbaled to reduce the data amount to be saved
        %fitswrite(single(phaseMap),[path_save,'ground_truth_phase_',idk,'_wfe_',sprintf('%.1f_', wfe),'nm.fits']);
        fitswrite(single(zCoefs),[path_save,'ground_truth_zernike_',idk,'_wfe_',sprintf('%.1f_', wfe),'nm.fits']);
        fitswrite(single(zCoefs_pyr),[path_save,'pyramid_based_zernike_',idk,'_wfe_',sprintf('%.1f_', wfe),'nm.fits']);
        fitswrite(single(pyr_frame),[path_save,'measurements_intensity_',idk,'_wfe_',sprintf('%.1f_', wfe),'nm.fits']);
        elapsed = toc(t_i);
        fprintf(['Time remaining :',num2str(elapsed*nAll/kIter/3600),' h\n']);
    end
end
    
tt = toc(t1);


