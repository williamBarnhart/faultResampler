% CALCSCALES  Calculates news smoothing scales for existing dislocations
%
%Usage
%Script, no function
%
%
% Edited June 24, 2010, by WDB
% Updated Nov 8, 2010 by WDB
% Cornell University
%
% Routine for calculating new smoothing scales for each dislocation using
% a Gaussian best-fit algorithm.  Inversion is first done using invertJRI
% for calculation of model resolution matrix, R.  
global covd2 numTri

nPatch  = length(patchstruct);

display(['Number of patches = ' num2str(nPatch)])


% Do inversion

switch rake_type
    case 'free'
        [gsmooth,G, Gg,mil]   = inversionFreeRake(patchstruct, resampstruct, Dnoise, lambdas, triId, flag);
                %[gsmooth,G, Gg,mil]   = invertJRI_freerake(patchstruct, resampstruct, Dnoise, lambdas, triId, 1);

        R       = Gg*G;
        r1      = R(1:nPatch,1:nPatch);
        r2      = R(nPatch+1:2*nPatch, nPatch+1:2*nPatch);
        R       = sqrt(r1.^2+r2.^2);

    case 'fixed'
        [gsmooth,G, Gg, slip, synth]   = inversionFixedRake(patchstruct, resampstruct, Dnoise, lambdas, triId, flag);
        R       = Gg*G;
        R       = R(1:nPatch, 1:nPatch);
end

Newscales=[];
tmpPatch=0;
for j=1:nfault
    xc      = mean([patchstruct(numTri(j)+1:numTri(j+1)).xfault]);
    yc      = mean([patchstruct(numTri(j)+1:numTri(j+1)).yfault]);
    zc      = mean([patchstruct(numTri(j)+1:numTri(j+1)).zfault]);
    
    %Calculate new length scale
    clear scales newscales
    for i=1:length(patchstruct(numTri(j)+1:numTri(j+1)))
        tmpPatch=tmpPatch+1;
        dist    =sqrt((xc-xc(i)).^2+(yc-yc(i)).^2+(zc-zc(i)).^2);
        scales1 =min(dist([1:i-1 i+1:length(patchstruct(numTri(j)+1:numTri(j+1)))]));
%         y2 =R(i+sum(numTri(1:j)),numTri(j)+1:numTri(j+1));
        y2 = R(tmpPatch,numTri(j)+1:numTri(j+1));
        b       = 0;
        x2      = [dist];
     
        id2     =find(y2>0);
        y2      =y2(id2);
        x2      =x2(id2);
        scales2= sqrt(sum(x2.^2.*y2)/sum(y2));
        scales(i)= max([scales1 scales2]);
        [out,resn,res]   = lsqnonlin('gausfun2',[scales(i)],[],[],options);
        newscales(i)=out(1);
        
    end
    Newscales=[newscales Newscales];
end
newscales=Newscales;

