function [green,greenX,greenY,greenZ,greenok]=make_green_meade_tri(patchstruct,datastruct,plotflag)
%ss +=left-lateral
%ds +=thrust

   
nu   = 0.25;
x    = [];
y    = [];
z    = [];
S    = [];

numfiles = length(datastruct);
Npatch   = length(patchstruct);

for i=1:numfiles;
  x    = [x; datastruct(i).X]';
  y    = [y; datastruct(i).Y]';
  S    = [S; [datastruct(i).S]'];
end
np = length(x);

green   = zeros(np,Npatch*2);
greenok = green;
if (nargout>1)
  greenX=green;
  greenY=green;
  greenZ=green;
end

%accept z values if given
if(isfield(datastruct,'z'))
    for i=1:numfiles;
        z    = [z; datastruct(i).Z];
    end
else
    z=x*0;
end

%h=waitbar(0,'Calculating Green''s Functions');
ts=0;
for j=1:2
 if(j==1)
     ds = 0;
     ss =-1;
 else
     ds = -1;
     ss = 0;
 end

 for i=1:Npatch
     id     = i+(j-1)*Npatch;
     vx=patchstruct(i).xfault;
     vy=patchstruct(i).yfault;
     vz=patchstruct(i).zfault;
     U        = CalcTriDisps(x',y',z',vx, vy, vz, nu, ss, ts, ds);
     green(:,id) = U.x.*S(:,1)+U.y.*S(:,2)-U.z.*S(:,3);
  
     if (nargout>1)
       greenX(:,id)= U.x;
       greenY(:,id)= U.y;
       greenZ(:,id)= -U.z;
   end
  
 end
end

if(plotflag)
    figure
    scatter(x,y,26,green(:,1),'filled')
    axis image,colorbar

end

