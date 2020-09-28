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

% Fitting and aliasing PSD
psdFit = fao.fittingPSD(fao.fx,fao.fy)/fao.atm.r0^(-5/3);

close all;
%% CREATING LINEAR HYPERBOLIC SAMPLING FOR INPUTS PARAMETERS - MAXIMIN - TRUNCATED 

r0_list = lhsdesign_modified(nr0,r0_min,r0_mean + r0_std);
v_list  = lhsdesign_modified(nv,v_min,v_mean + v_std);
n_list  = lhsdesign_modified(nn,n_min,n_mean + n_std);
nAll    = nr0*nv*nn

%% LOOP
kIter = 0;
elapsed = 0;

for iv = 1:1
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
            % update the noise value
            fao.noiseVariance = n_list(in);
                        
            % Get the residual phase PSD
            psdAO = psdFit*atm.r0^(-5/3);
            psdAO(index) = psdAO(index) + reshape(fao.noisePSD(fao.fx,fao.fy) + psdServo,[],1);
            
            % Get the phase map
            phaseMap = tools.crop( pupZeroPad .* real( fft2( fftshift( sqrt(psdAO).*randn(size(psdAO)) ) ) ),nPxPup);
            
            % propagate through the pyramid
            ngs.phase = 0;
            ngs.phase = phaseMap;
            ngs = ngs*tel*pyr;
            
            % record outputs
            fitswrite(single(phaseMap),[path_save,'ground_truth_phase_r0_',num2str(r0_list(ir0)),'_v_',num2str(v_list(iv)),'_noise_',num2str(n_list(in)),'.fits']);
            fitswrite(single(pyr.camera.frame),[path_save,'measurements_intensity_r0_',num2str(r0_list(ir0)),'_v_',num2str(v_list(iv)),'_noise_',num2str(n_list(in)),'.fits']);
            elapsed = elapsed + toc(t_i);
            fprintf(['Time remaining :',num2str(elapsed/kIter*(nAll-kIter)/3600),' h\n']);
            
        end
    end
end
close all;