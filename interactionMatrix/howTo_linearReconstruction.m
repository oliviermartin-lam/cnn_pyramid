% ---------------------------------------------
%   How to: Linear Reconstruction with Pyramid
%
%   Pyramid - Deep Learning Project
% ---------------------------------------------

clear all
close all
clc
wfsModel = 'diff';

%% Definition of the atmosphere 
r0 = 1;
atm = atmosphere(photometry.V,r0,30,...
    'altitude',[0,4,10]*1e3,...
    'fractionnalR0',[0.7,0.25,0.05],...
    'windSpeed',[5,10,20],...
    'windDirection',[0,pi/4,pi]);

%% Definition of the telescope
nPx = 60;
tel = telescope(8,...
    'fieldOfViewInArcMin',2.5,...
    'resolution',nPx,...
    'samplingTime',1/500);

%% Zernike Base
nModes = 50;
zer = zernike(2:nModes+1,tel.D,'resolution',nPx);
modes = zer.modes;

%% Definition of a calibration source
ngs = source('wavelength',photometry.V);

%% Definition of the wavefront sensor
nLenslet = nPx;
modu = 3;
wfs = pyramid(nLenslet,nPx,'c',2,'modulation',modu);
% ----- NOISE --------------
ngs = ngs.*tel*wfs;
wfs.INIT
+wfs;

%% Interaction Matrix - FULL FRAME
iMat = interaction_matrix(ngs,tel,wfs,modes);
commandMatrix = pinv(iMat);

%% Reconstruction of a given phase

% Make Phase - Here with turbulence
tel = tel+atm;
ngs = ngs.*tel;
phase = ngs.phase;

% Decompostion on Zernike
phase2zernike = pinv(modes);% Zernike reconstructor
zernike_phase = phase2zernike*phase(:);%Decompostion on zernike

% ----- Propagation on Pyramid -----
tel = tel - atm;% remove atmosphere
ngs = ngs.*tel;
% Flat for FULL FRAME
ngs = ngs*wfs;
I_0 = wfs.camera.frame./sum(wfs.camera.frame(:));
% Phase through Pyramid
ngs = ngs.*tel;
ngs.phase = phase;% Put phase to estimate
ngs = ngs*wfs;%propagation on Pyramid
fullFrame = wfs.camera.frame./sum(wfs.camera.frame(:))-I_0;

% Reconstruction
zernike_phase_estimated = commandMatrix*fullFrame(:);

%% PLOT
figure;
plot(zernike_phase);
hold on;
plot(zernike_phase_estimated);
grid
legend('Phase','Estimated Phase - linear reconstructor');





