global b x2 y2 covd2

% FAULTRESAMPLER    Tool for generating optimally discretized model fault
% planes
%
% Usage
%  Script, no function
%  
%  Edited June 24, 2010
%  Cornell University
%
%  Updates
%   November 9, 2010: Added lines to calcScales.m to allow for free rake
%                       inversion.  Also added invertJRI_freerake.m for 
%                       free rake inversions.  Changed citation to
%                       appropriate journal citation.
%  
% This program generates fault models for inversions of slip from GPS and
% InSAR observed surface displacements.  It then conducts a final inversion
% of slip using the optimal fault model. 
%
% Citation:  Barnhart, W.D., R.B. Lohman (2010) Automated fault model
% discretization for inversions for coseismic slip distributions. Journ. 
% Geophys. Res. V.115, B10419
%

% Turn warnings off singular matrices and delaunay triangulations as first
% iterations will have singular Green's matrices or too few model
% parameters.
% warning off MATLAB:nearlySingularMatrix
% warning off MATLAB:delaunayn:DuplicateDataPoints
% warning('off','Duplicate data points have been detected and removed. The Triangulation indices are defined with respect to the unique set of points in DelaunayTri property X.')
% 
warning off

loadResampData


makeStartFault_multi

OK   = 0;
iter= 1;

while OK<1
    calcScales_multi
    makeNewMesh_multi
    checkScales_multi
end

% Do inversion of final 
nPatch  = length(patchstruct);


switch rake_type
    case 'free'
        [gsmooth,G, Gg, slip, synth, mil]   = inversionFreeRake(patchstruct, resampstruct, Dnoise, lambdas, triId, 0);
%                 [gsmooth,G, Gg, slip, synth, mil]   = invertJRI_freerake(patchstruct, resampstruct, Dnoise, lambdas, triId, 0);

    case 'fixed'
        [gsmooth,G, Gg, slip, synth]   = inversionFixedRake(patchstruct, resampstruct, Dnoise, lambdas, triId, 0);
end

rotateFinal

figure
patch([patch_new1.yfault], [patch_new1.zfault], slip(1:length(patch_new))')
axis image
set(gca,'ydir', 'reverse')
colorbar
title('Final Slip Distribution')

figure
subplot(3,1,1)
patch(boxx,boxy, data');
axis image
shading flat
colorbar
title('Data')

subplot(3,1,2)
patch(boxx,boxy, synth')
axis image
shading flat
colorbar
title('Synthetic')

subplot(3,1,3)
patch(boxx,boxy, data'-synth')
axis image
shading flat
colorbar
title('Data Misfit')


    
