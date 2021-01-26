%% ATMOSPHERE
L0          = 25;               % outer scale in meters
r0_min      = 0.01;r0_max=0.11; % min and max values of r0 in meters
v_min       = 3;v_max  = 30;    % min and max values of wind speed in m/s

%% GUIDE STAR
photoNGS    = photometry.I1; %700nm
magNGS      = 5;
%% TELESCOPE
D           = 1.5;              % telescope diameter in meter
cobs        = 0;                % telescope central obstruction in percent
nPxPup      = 64;               % number of pixels to describe the pupil

%% PYRAMID PARAMETERS

nLenslet    = 16;               % number of lenslets for the pyramid
pyrBinning  = 1;                % binning rate on the pyramid detector
pyrMod      = 3;                % Modulation rate on top of the pyramid in lambda/D units
Samp        = 2;                % OVer-sampling factor
resAO       = 2*nLenslet+1;     % number of pixels to describe the AO-correction area in the PSD
fovInPixel  = nPxPup*2*Samp;    % number of pixel to describe the PSD
nTimes      = fovInPixel/resAO;
ron         = 0.1;              % read-out noise in e-
QE          = 0.6 * 0.8;        % total optical throughput * Quantum efficiency

%% AO LOOP USED FOR THE FOURIER GENERATOR
loopGain    = 0.5;              % AO loop main gain
sampTime    = 3.5e-3;           % AO sampling time in seconds
latency     = 2.5*sampTime;     % AO loop latency in seconds;
n_min       = 0.1;n_max= 1.5;   % noise level in rad^2

%% GROUND TRUTH
nAll        = 1250000;          % Number of simulations 1250000
genType     = 'Zernike';        %can be Zernike, or Fourier
jIndex      = [4];          % Noll's indexes of the Zernike modes to be simulated for the Zernike Generator - 1 is piston
nZern       = 50;               % number of Zernike to be reconstructed for the Fourier generator
zStdMax     = [500];            % Vector of size (1,numel(jIndex)) containing Std of the Zernike amplitude in nm
zMean       = [0];              % Vector of size (1,numel(jIndex)) containing mean of the Zernike amplitude in nm
zDistrib    = 'Normal';         % Normal -> N(zMean, zStdMax) or uniform -> [zMean-5*zStdMax,zMean+5*zStdMax]
