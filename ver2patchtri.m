function patchstruct =   ver2patchtri(faultstruct, p, t)

% VER2PATCHTRI  Make a patch structure usable in Green's funtion algorithms
% for triangular dislocations
%
% Usage
% [patchstruct] = ver2patchtri(faultstruct, p, t)
%
% Variables:
%   FAULTSTRUCT - geometric information about a single planar fault segment
%
%   P - X-Y gird coordinates of triangle dislocation vertices.  Same as
%   triCoords in CHECKSCALES and MAKENEWMESH
%
%   T - Index terms for triangles with coordinates defined in P.  Same as
%   triId in CHECKSCALES and MAKENEWMESH
%
%   Citation:Barnhart, W. D., and R. B. Lohman (2010), Automated fault 
%     model discretization for inversions for coseismic slip distributions, 
%     J. Geophys. Res., 115, B10419, doi:10.1029/2010JB007545.




numvert = length(p);
nt      = size(t,1);
L       = faultstruct.L;
W       = faultstruct.W;
strike  = faultstruct.strike;
dip     = faultstruct.dip;
zt      = faultstruct.zt;
vertices= faultstruct.vertices;
xc      = mean(vertices(1,:));
yc      = mean(vertices(2,:));
x       = zeros(numvert, 1); %line of strike
y       = p(:,1); %Normal to strike
z       = p(:,2); %Depth Dip direction
p1      = [x y z];

%Strike rotation CW around Z axis
T2      = [-sind(dip) 0 cosd(dip); 0 1 0; cosd(dip) 0 sind(dip)];
T       = [cosd(strike) -sind(strike) 0; sind(strike) cosd(strike) 0; 0 0 1];
p2      = p1*T2;

%Dip rotation CW around line of strike
p3      = p2*T;
p3(:,1) = p3(:,1)+xc;
p3(:,2) = p3(:,2)+yc;
p3(:,3) = p3(:,3)+zt;

% plot3(p3(:,1), p3(:,2), p3(:,3), 'k*')

for i=1:nt
    patchstruct(i).xfault= [p3(t(i,:), 1)];
    patchstruct(i).yfault= [p3(t(i,:), 2)];
    patchstruct(i).zfault= [p3(t(i,:), 3)];
end
% patchstruct.triId     = t;
% patchstruct.triCoords = p;
