function smooth= triSmooth(t)

% TRISMOOTH   Laplacian smoothing matrix for triangles
%
% Usage
% [smooth] = triSmooth(t)
%
% Edited June 24, 2010 by WDB
% Cornell Univeristy
%
% Routine generates a Laplacian smoothing matrix for higher level Tikhonov
% regularization.  Smoothing of triangles is not weighted by their area
% since it is assumed that triangle size gradients vary smoothly.
%
% Variables
% T - Index terms for triangles with coordinates defined in P.  Same as
% triId in CHECKSCALES and MAKENEWMESH


smooth = eye(length(t));
for i=1:length(t)
    clear common
    common= [];
    a   = t(i,:);
    for j= 1:length(t)
        b=t(j,:);
        tmp = intersect(a,b);
        if length(tmp)==2
           common= [common j];
        end
    end
    nneighbor= length(common);
    
    if nneighbor==1
        id= common;
        smooth(i,id)= -1;
    elseif nneighbor==2
        smooth(i,common)= -1/2;
    else
        smooth(i,common)=-1/3;
    end

end

            
    