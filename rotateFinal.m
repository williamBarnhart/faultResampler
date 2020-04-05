% rotateFinal
clear patch_new patch_new1
    strike=faultstruct.strike;
    dip=faultstruct.dip;


rot1= [cosd(-strike) sind(-strike) 0; -sind(-strike), cosd(-strike) 0; 0 0 1];

for i=1:length(patchstruct)
    
    % clear tmp
    tmp= rot1*[[patchstruct(i).xfault]';[patchstruct(i).yfault]';[patchstruct(i).zfault]'];
    % tmp2 = rot2*tmp;
    % tmp2= rot2*tmp;
    patch_new1(i).xfault=tmp(1,:)';
    patch_new1(i).yfault=tmp(2,:)';
    patch_new1(i).zfault=tmp(3,:)';
end


%Strike rotation

rot2= [cosd(dip-90) 0 sind(dip-90); 0 1 0; -sind(dip-90) 0 cosd(dip-90)];
for i=1:length(patchstruct)
    tmp2=rot2*[[patch_new1(i).xfault]';[patch_new1(i).yfault]';[patch_new1(i).zfault]'];
    patch_new(i).xfault=tmp2(1,:)';
    patch_new(i).yfault=tmp2(2,:)';
    patch_new(i).zfault=tmp2(3,:)';
end
