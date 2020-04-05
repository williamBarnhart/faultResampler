function [gsmooth, G, Gg, mil, synth] = inversionFixedRake(patchstruct, resampstruct,Dnoise, lambdas, triId, flag)

% INVERTJRI  Linear inversion using jRi
%
% Usage
% [gsmooth, G, Gg, mil, synth] = invertJRI(patchstruct, resampstruct,Dnoise, lambdas, triId, flag)
%
%
%   Edited June 24, 2010, by WDB
%   Cornell University
%
%   Citation:
%   Barnhart, W. D., and R. B. Lohman (2010), Automated fault model
%   discretization for inversions for coseismic slip distributions, J.
%   Geophys. Res., 115, B10419, doi:10.1029/2010JB007545.
%
%
%Algorithm to invert for distributed slip using Laplacian smoothing and
% Tikhonov regularization.  Smoothing parameter is chosen using the approximate jRi
% method
%
% Variables:
% PATCHSTRUCT - data structure with information about the fault
% orientation and geometry of individual dislocations
%
% RESAMPSTRUCT- data structure with information about the surface
% data locations, values, and LOS vector
%
% DNOISE - data matrix weighted by the inverse of the cholesky
% factorization
%
% ALPHAS - potential smoothing factors for jRi to choose from
%
% TRIID -
%
% FLAG - plotting option.  FLAG=0 no figs, FLAG=1 makes figs

global rampg Cdinv nramp rake covd2 smooth_method reg_method data_type faultstruct

options     = optimset('LargeScale', 'off', 'MaxIter', 1000);
nPatch      = length(patchstruct);
np          = length([resampstruct.data]);

[green]     = make_green_meade_tri(patchstruct, resampstruct, 0);
g1          = green(:,1:nPatch)';               % SS component
g2          = green(:,nPatch+[1:nPatch])';      % DS component
green       = cosd(rake)*g1-sind(rake)*g2;
% green       = cosd(rake)*g1+sind(rake)*g2;    %same convention as NA

green       = green;
G           = Cdinv*[green; rampg]';             % Weighted Greens Functions
D           = [Dnoise; zeros(nPatch,1)];

switch smooth_method
    case 'laplacian'
        smooth      = triSmooth(triId);                 %Laplacian Smoother
    case 'mm'
        smooth      = eye(nPatch);                      %Minimum moment
end

A           = blkdiag(-1*eye(nPatch),zeros(nramp));
B           = zeros(nPatch+nramp,1);



%%%%% jRi
for i=length(lambdas):-1:1
    mil         = zeros(1,nPatch+nramp);
    
    lambda              = lambdas(i);
    gsmooth             = [G; lambda*smooth zeros(nPatch, nramp)];
    [mil, resnorm, ril] = lsqlin(gsmooth,D,A,B,[],[],[],[],mil, options);
    Gg                  = inv(gsmooth'*gsmooth)*G';
    N                   = G*Gg;
    M                   = [eye(np) -N];
    iRi                 = sum(ril(1:np).^2)/np;
    covresjRi           = M*M';
    covresiRi           = M*covd2*M';
    jRin                = mean(diag(covresjRi));
    iRin                = mean(diag(covresiRi));
    oRo_approx          = iRi-iRin;
    jRi(i)              = oRo_approx+jRin;
    r_norm(i)           = iRi;
    m_norm(i)           = (smooth*mil(1:nPatch))'*(smooth*mil(1:nPatch));
end


switch reg_method
    case 'jRi'
        id              = find(jRi==min(jRi));
    case 'lcurve'
        h   = figure;
        subplot(1,2,1)
        plot(sqrt(r_norm),sqrt(m_norm),'-o')
        title('Choose L-Curve corner ID')
        xlabel('Data Norm')
        ylabel('Model Norm')
        
        subplot(1,2,2)
        semilogx(lambdas, jRi, '.-')
        
        id  = input([' \n'...
            '\n'...
            'Choose ID number of corner point, starting from the right \n','s']);
        close(h)
end


lambda              = lambdas(id);
gsmooth             = [G; lambda*smooth zeros(nPatch, nramp)];
Gg                  = inv(gsmooth'*gsmooth)*G';
[mil,resnorm,ril]   = lsqlin(gsmooth, D, A, B, [],[],[],[],mil,options);
synth               = green'*mil(1:nPatch);
[m0,mw]             = calcMoment(patchstruct, mil(1:nPatch), 'tri')
slip                = mil(1:nPatch);

%%%% Plot Stuff
if flag==0
    return
else
    switch data_type
        case 'InSAR'
            figure
            
            subplot(3,2,1)
            semilogx(lambdas, jRi, '.-')
            hold on
            plot(lambdas(id), jRi(id), 'ro')
            axis tight
            xlabel('lambda')
            ylabel('jRi')
            title('jRi Curve')
            
            subplot(3,2,2)
            plot(sqrt(r_norm), sqrt(m_norm), '.-')
            hold on
            plot(sqrt(r_norm(id)), sqrt(m_norm(id)), 'ro')
            axis tight
            title('L-curve')
            
            rotateFinal
            subplot(3,2,3)
            patch([patch_new1.yfault], [patch_new1.zfault], slip')
            axis image
            colorbar
            set(gca,'ydir','reverse')
            title('Inverted Slip')
            
            subplot(3,2,4)
            scatter([resampstruct.X], [resampstruct.Y], 20, [resampstruct.data]', 'filled')
            axis image
            colorbar
            title('Data')
            
            subplot(3,2,5)
            scatter([resampstruct.X], [resampstruct.Y], 20, synth', 'filled')
            axis image
            colorbar
            title('Model')
            
            subplot(3,2,6)
            scatter([resampstruct.X], [resampstruct.Y], 20, [resampstruct.data]'-synth', 'filled')
            axis image
            colorbar
            title('Misfit')
        case 'GPS'
            %             odds = [1:2:np];
            %             evens= [2:2:np];
            S = [resampstruct.S];
            odds = find(S(1,:)==1);
            evens= find(S(2,:)==1);
            
            figure
            subplot(3,2,1)
            semilogx(lambdas, jRi, '.-')
            hold on
            plot(lambdas(id), jRi(id), 'ro')
            axis tight
            xlabel('lambda')
            ylabel('jRi')
            title('jRi Curve')
            
            subplot(3,2,2)
            plot(sqrt(r_norm), sqrt(m_norm), '.-')
            hold on
            plot(sqrt(r_norm(id)), sqrt(m_norm(id)), 'ro')
            axis tight
            title('L-curve')
            
            rotateFinal
            subplot(3,2,3)
            patch([patch_new.yfault], [patch_new.zfault], slip')
            colorbar
            set(gca,'ydir','reverse')
            title('Inverted Slip')
            
            subplot(3,2,4)
            quiver([resampstruct.X(odds)],[resampstruct.Y(odds)],[resampstruct.data(odds)],[resampstruct.data(evens)], 1);
            axis image
            title('Data');
            
            subplot(3,2,5)
            quiver([resampstruct.X(odds)],[resampstruct.Y(odds)],[resampstruct.data(odds)],[resampstruct.data(evens)], 1);
            hold on
            quiver([resampstruct.X(odds)],[resampstruct.Y(odds)],synth(odds),synth(evens), 1,'r');
            axis image
            title('Model')
            
            subplot(3,2,6)
            quiver([resampstruct.X(odds)],[resampstruct.Y(odds)],[resampstruct.data(odds)]-synth(odds),[resampstruct.data(evens)]-synth(evens), 1);
            axis image
            title('Misfit')
            figure
            
        case 'Mixed'
            
            subplot(3,2,1)
            semilogx(lambdas, jRi, '.-')
            hold on
            plot(lambdas(id), jRi(id), 'ro')
            axis tight
            xlabel('lambda')
            ylabel('jRi')
            title('jRi Curve')
            
            subplot(3,2,2)
            plot(sqrt(r_norm), sqrt(m_norm), '.-')
            hold on
            plot(sqrt(r_norm(id)), sqrt(m_norm(id)), 'ro')
            axis tight
            title('L-curve')
            
            rotateFinal
            subplot(3,2,3)
            patch([patch_new.yfault], [patch_new.zfault], slip')
            colorbar
            set(gca,'ydir','reverse')
            title('Inverted Slip')
            
            subplot(3,2,4)
            scatter([resampstruct.X], [resampstruct.Y], 20, [resampstruct.data]', 'filled')
            axis image
            colorbar
            title('Data')
            
            subplot(3,2,5)
            scatter([resampstruct.X], [resampstruct.Y], 20, synth', 'filled')
            axis image
            colorbar
            title('Model')
            
            subplot(3,2,6)
            scatter([resampstruct.X], [resampstruct.Y], 20, [resampstruct.data]'-synth', 'filled')
            axis image
            colorbar
            title('Misfit')
    end
end