% CHECKSCALES   Check termination criterion given the predefined tolerance
% limit
%
% Usage:
% Script, no function
%
% Editd June 24, 2010
% Cornell University
%
% This routine checks the number of dislocations in the newly formed mesh,
% saved as patchstructnew.  If the percent difference between the number of
% dislocations in the current mesh and the previous iteration are less than
% the tolerance, resampling algorithm is terminated and the most recent
% discretization is made the final discretization

nTri    = size(triId, 1);

% clear patchstructnew
% patchstructnew = ver2patchtri(faultstruct, triCoords, triId);
npatch  = length(patchstructnew);
nDiff   = npatch*tol;

if (abs(nPatch-npatch)<=nDiff)
    OK=1;
else
    OK=0;
end

patchstruct = patchstructnew;
display(['Number of patches = ' num2str(length(patchstruct))])
