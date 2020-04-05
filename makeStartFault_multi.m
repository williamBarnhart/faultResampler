global nfault

nfault=length(faultstruct)
patchstruct=[];
triId= [];
numTri=[];

for i=1:nfault
    Faultstruct=faultstruct(i);
   L=Faultstruct.L;
   dip=Faultstruct.dip;
   W=Faultstruct.W;
   zt=Faultstruct.zt;


node            = [-L/2 0; L/2 0; L/2 W; -L/2 W];
yc              = rand(100,1)*L-L/2;
zc              = rand(100,1)*W;
newscales       = 3000*ones(size(yc)); %Adjust constant for different sized faults
hdata.fun       = @hfun1;
hdata.args      = {newscales, yc, [(zc-zt)./sind(dip)]};
[triCoords, TriId]= mesh2d(node,[], hdata, options3);

nTriId(i)= length(TriId);
sums=sum(nTriId(1:i));
numTri= [numTri sums];

close
Patchstruct     = ver2patchtri(Faultstruct, triCoords, TriId);
patchstruct=[patchstruct Patchstruct];
triId= [triId; TriId];
end
numTri=[0 numTri];


display(['Number of patches = ' num2str(length(patchstruct))])

