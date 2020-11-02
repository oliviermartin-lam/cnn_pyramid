function [ iMat ] = interaction_matrix(ngs,tel,wfs,modalBasis)
% Compute intereaction matrix with push pull method around flat wavefront
% - signal: full frame
    amplitude = 0.1;%Need to be small (for noise-free esystem)
    modalBasis = reshape(modalBasis,tel.resolution,tel.resolution,size(modalBasis,2));%Reshape modes
    % ------- PUSH ---------
    ngs = ngs.*tel;
    ngs.phase = amplitude*modalBasis;
    ngs = ngs*wfs;
    camera = reshape(wfs.camera.frame,2*wfs.c*wfs.nLenslet,2*wfs.c*wfs.nLenslet,[]);
    % Normalisation
    sp= camera./sum(sum(camera(:,:,1)));
    
    % ------- PULL ---------
    ngs = ngs.*tel;
    ngs.phase = -amplitude*modalBasis;
    ngs = ngs*wfs;     
    camera = reshape(wfs.camera.frame,2*wfs.c*wfs.nLenslet,2*wfs.c*wfs.nLenslet,[]);
    % Normalisation
    sm= camera./sum(sum(camera(:,:,1)));
    
    % PUSH-PULL COMPUTATION
    iMat=0.5*(sp-sm)/amplitude;
    
    % Reshape in 2-D
    iMat = reshape(iMat,(2*wfs.c*wfs.nLenslet)^2,[]);
  
end