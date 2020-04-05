function savestructstruct = writeDatastruct(X,Y,data,S,cov,zone,type);
%%% WRITEDATASTRUCT.M
%%% Write data points to a Matlab structure that can be read by
%%% loadResampData.m
%%%
%%% X,Y: X and Y locations of a data point, in UTM. Variables should be
%%%     vectors of size [npx1]. np is the total number of data points
%%%
%%% data: Magnitude of displacement at point (X,Y). If using InSAR, it should be
%%%     [np x 1]. For GPS, it should either be [np x 2] (E-W, N-S), or [np x 3]
%%%     (E-W, N-S, vertical);
%%%
%%% S: Look vector. If using InSAR, S should be [3 x np]. If using GPS,
%%% enter [];
%%%
%%% cov: Data covariance matrix. For InSAR, this should be sized [np x np]
%%%     and may either be a diagonal matrix (variances on the diagonal), or a
%%%     full covariance matrix. For GPS, this should be a matrix that is
%%%     either [np x 2] or [np x 3], with East errors in the first column,
%%%     north errors in the second, and vertical (if applicable) in the
%%%     third
%%%
%%% zone: UTM zone identifier, e.g. '10S'
%%%
%%% type: Either 'InSAR' or 'GPS'

switch type
    case 'InSAR'
        np      = length(data);
        fill    = zeros(1,np);      % This is a dumby variable
        for     i = 1:np
            savestruct.data(i).X        = X(i);
            savestruct.data(i).Y        = Y(i);
            savestruct.data(i).S        = S(:,i);
            savestruct.data(i).data     = data(i);
            savestruct.covstruct.cov    = cov;
            %fill in dumby variables
            savestruct.data(i).boxx     = zeros(5,1);
            savestruct.data(i).boxy     = zeros(5,1);
        end
        savestruct.zone = zone;
        savestruct.np   = np;
        savestruct.type = type;
        
    case 'GPS'
        [np, nComponents] = size(data);
        if nComponents ==1
            error('Only one component given, need 1 or 2 more components');
        elseif nComponents == 2
            S       = zeros(3,np);
            dE      = data(:,1);
            dN      = data(:,2);
            eE      = cov(:,1);
            eN      = cov(:,2);
            for i = 1:np
                X(2*i-1,1)      = X(i);
                Y(2*i-1,1)      = Y(i);
                X(2*i,1)        = X(i);
                Y(2*i,1)        = Y(i);
                S(1, i*2-1)     = 1;
                S(2, i*2)       = 1;
                data(2*i-1,1)   = dE(i);
                data(2*i,1)     = dN(i);
                covd(2*i-1,1)   = eE(i);
                covd(2*i,1)     = eN(i);
            end
            
            covd    = diag(covd);
            np      =np*2;
            fill    = zeros(1,np);
            for     i = 1:np
                savestruct.data(i).X        = X(i);
                savestruct.data(i).Y        = Y(i);
                savestruct.data(i).S        = S(:,i);
                savestruct.data(i).data     = data(i);
                savestruct.covstruct.cov    = covd;
                %fill in dumby variables
                savestruct.data(i).boxx     = zeros(5,1);
                savestruct.data(i).boxy     = zeros(5,1);
            end
            savestruct.zone = zone;
            savestruct.np   = np;
            savestruct.type = type;
            
        else
            S   = zeros(3,np);
            dE  = data(:,1);
            dN  = data(:,2);
            dZ  = data(:,3);
            eE  = cov(:,1);
            eN  = cov(:,2);
            eZ  = cov(:,3);
            for i = 1:np
                X(3*i-2,1)      = X(i);
                Y(3*i-2,1)      = Y(i);
                X(3*i-2,1)      = X(i);
                Y(3*i-1,1)      = Y(i);
                X(3*i,1)        = X(i);
                Y(3*i,1)        = Y(i);
                S(1, i*3-2)     = 1;
                S(2, i*3-1)     = 1;
                S(3, i*3)       = 1;
                data(3*i-2,1)   = dE(i);
                data(3*i-1,1)   = dN(i);
                data(3*i,1)     = dZ(i);
                covd(3*i-2,1)   = eE(i);
                covd(3*i-1,1)   = eN(i);
                covd(3*i,1)     = eZ(i);
            end
            
            covd    = diag(covd);
            np      =np*3;
            fill    = zeros(1,np);
            for     i = 1:np
                savestruct.data(i).X        = X(i);
                savestruct.data(i).Y        = Y(i);
                savestruct.data(i).S        = S(:,i);
                savestruct.data(i).data     = data(i);
                savestruct.covstruct.cov    = covd;
                %fill in dumby variables
                savestruct.data(i).boxx     = zeros(5,1);
                savestruct.data(i).boxy     = zeros(5,1);
            end
            savestruct.zone = zone;
            savestruct.np   = np;
            savestruct.type = type;
        end
end
