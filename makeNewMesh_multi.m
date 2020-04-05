% MAKENEWMESH   Use new scales to calculate a new triangular mesh using
% MESH2D
%
% Usage:
% Script, no function
% 
% Edited June 24, 2010 by WDB
% Cornell University
numtri=[];
clear triCoords triId
triId=[];
patchstructnew=[];
for i=1:nfault
    
    
    node = [-faultstruct(i).L/2 0; faultstruct(i).L/2 0; faultstruct(i).L/2 faultstruct(i).W; -faultstruct(i).L/2 faultstruct(i).W];
x1          = faultstruct(i).vertices(1,1);
y1          = faultstruct(i).vertices(2,1);
dl          = sqrt((xc-x1).^2+(yc-y1).^2)-L/2;
hdata.args  = {newscales(numTri(j)+1:numTri(j+1)), dl, [(zc-[faultstruct(i).zt])./sind(faultstruct(i).dip)]};



[triCoords, TriId] = mesh2d(node, [], hdata, options3);
close

ntriid(i)=length(TriId);
sums= sum(ntriid(1:i));
numtri=[numtri sums];
triId= [triId; TriId];
Patchstruct     = ver2patchtri(faultstruct(i), triCoords, TriId);
patchstructnew= [patchstructnew Patchstruct];
end
numTri=[0 numtri];


    nTri    = length(patchstructnew);

% clear patchstructnew
% Patchstructnew = ver2patchtri(faultstruct, triCoords, triId);






