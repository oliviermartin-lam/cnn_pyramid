
%% GUIDE STAR
photoNGS    = photometry.V0;

%% TELESCOPE
D           = 1;  % telescope diameter in meter
cobs        = 0;  % telescope central obstruction in percent
nPxPup      = 60; % number of pixels to describe the pupil

%% PYRAMID PARAMETERS

nLenslet    = nPxPup; % number of lenslets for the pyramid
pyrBinning  = 1;  % binning rate on the pyramid detector
pyrMod      = 3;  % Modlulation rate on top of the pyramid in lambda/D units
noiseVar    = 1;  % Noise variance in rd^2
Samp        = 2;  % OVer-sampling factor
resAO       = nPxPup/2 + 1; % number of pixels to describe the AO-correction area in the PSD
fovInPixel  = nPxPup*2*Samp; % number of pixel to describe the PSD
nTimes      = fovInPixel/resAO;


%% AO LOOP
loopGain    = 0.5;  % AO loop main gain
sampTime    = 3.3e-3; % AO sampling time in seconds
latency     = 1*sampTime; % AO loop latency in seconds;

%% ATMOSPHERE
L0      = 25;

nr0     = 50;
r0_min  = 0.01;
r0_mean = 0.08;
r0_std  = 0.03;

nv      = 50;
v_min   = 1;
v_mean  = 10;
v_std   = 10;

nn      = 50;
n_min   = 0;
n_mean  = 0.1;
n_std   = 1;