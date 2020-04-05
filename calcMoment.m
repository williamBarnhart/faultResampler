function [M0, Mw]= calcMoment(patchstruct, slip, type)

slip= abs(slip);
switch type
    case 'quad'
        id          = find(slip>0);
        for j=1:length((id))
            
         area        = sum([patchstruct(id(j)).L].*[patchstruct(id(j)).W])*1e4; %area in cm^2
         gamma(j)       = area*sum(slip(id(j))*1e2);
        end
         eta         = sum(gamma);
    case 'tri'
        id          = find(slip>0);
        if id>0
        for i= 1:length(id)
            j= id(i);
            first(i)       = (det([[patchstruct(j).xfault(1)] [patchstruct(j).xfault(2)] [patchstruct(j).xfault(3)]; [patchstruct(j).yfault(1)] [patchstruct(j).yfault(2)] [patchstruct(j).yfault(3)]; 1 1 1]));
            second(i)      = (det([[patchstruct(j).yfault(1)] [patchstruct(j).yfault(2)] [patchstruct(j).yfault(3)]; [patchstruct(j).zfault(1)] [patchstruct(j).zfault(2)] [patchstruct(j).zfault(3)]; 1 1 1]));
            third(i)       = (det([[patchstruct(j).zfault(1)] [patchstruct(j).zfault(2)] [patchstruct(j).zfault(3)]; [patchstruct(j).xfault(1)] [patchstruct(j).xfault(2)] [patchstruct(j).xfault(3)]; 1 1 1]));
            gamma(i)       = (1/2*sqrt(first(i)^2+second(i)^2+third(i)^2)*1e4)*slip(j)*1e2;
%                         gamma(i)       = (1/2*sqrt(first(i)^2+second(i)^2+third(i)^2)*1e4)*slip(i)*1e2;

        end
        
        eta     = sum(gamma);
        
        else
        eta=0;
        end
end

%totalslip   = sum(slip(id))*1e2; %Slip in cm, scale if necessary
% mu          = 3.0e11; %dyne/cm^2
mu =4.0e11;
M0          = eta*mu;
Mw          = log10(M0)/1.5-10.73;