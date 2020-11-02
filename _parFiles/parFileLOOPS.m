%% ATMOSPHERE
L0      = 25;
r0_min  = 0.01;r0_max  = 0.11;
v_min   = 3;v_max  = 30;

%% GUIDE STAR
photoNGS    = photometry.I1;

%% TELESCOPE
D           = 1.5;  % telescope diameter in meter
cobs        = 0;  % telescope central obstruction in percent
nPxPup      = 64; % number of pixels to describe the pupil

%% PYRAMID PARAMETERS

nLenslet    = 16; % number of lenslets for the pyramid
pyrBinning  = 1;  % binning rate on the pyramid detector
pyrMod      = 3;  % Modlulation rate on top of the pyramid in lambda/D units
Samp        = 2;  % OVer-sampling factor
resAO       = 2*nLenslet+1; % number of pixels to describe the AO-correction area in the PSD
fovInPixel  = nPxPup*2*Samp; % number of pixel to describe the PSD
nTimes      = fovInPixel/resAO;

%% AO LOOP
loopGain    = 0.5;  % AO loop main gain
sampTime    = 3.5e-3; % AO sampling time in seconds
latency     = 2.5*sampTime; % AO loop latency in seconds;
n_min       = 0.1;n_max = 1.5;
nZern       = 50;

