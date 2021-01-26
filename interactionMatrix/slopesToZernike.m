function S2Z = slopesToZernike(tel,wfs,sref,varargin)
inputs = inputParser;
inputs.addRequired('tel', @(x) isa(x,'telescope'));
inputs.addRequired('wfs', @(x) isa(x,'shackHartmann') || isa(x, 'pyramid'));
inputs.addRequired('sref',@(x) isa(x,'source'));
inputs.addParameter('nModes', 200,@isnumeric);
inputs.addParameter('nThresholded',0,@isnumeric);
inputs.addParameter('amp',sref.wavelength/40,@isnumeric);
inputs.parse(tel,wfs,sref,varargin{:});
nModes          = inputs.Results.nModes;
nThresholded    = inputs.Results.nThresholded;
amp             = inputs.Results.amp;

%1\ Init
wfs.camera.readOutNoise = 0;
wfs.camera.photonNoise = false;
%wfs.slopes = 0*wfs.slopes;

%2\ Define Zernike modes, piston-excluded
zer   = zernike(2:nModes+1,'resolution',tel.resolution,'D',tel.D);

%3\ Calibration
dm                      = deformableMirror(nModes,'modes',zer,'resolution',tel.resolution);
sref                    = sref.*tools.duplicateTelescope(tel); %duplicate the telescope so as to avoid including the atmosphere durinf the calibration process
dmCalib                 = calibration(dm,wfs,sref,amp,1);
dmCalib.nThresholded    = nThresholded;

% Zernike reconstruction matrix in pixels/nm
S2Z = -2e9*dmCalib.M;
