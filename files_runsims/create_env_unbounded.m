function [] = create_env_unbounded(sx,sy,nbins_env,b)
% written by Allison K. Shaw (ashaw@umn.edu)
% updated: June 2018
%
% USAGE: [] = create_env_unbounded(sx,sy,nbins_env,b)
%
% Create the environment to use as a backdrop for the dispersal IBM.
% Here, the environment is a square of viable patches set of patches with
% wrapping boundaries.
% Only need to run this once - it creates a mat file that can be uploaded
%   later within IBM code.
%
% INPUTS:
%   sx = number of sites in the x-dimension of the environment
%   sy = number of sites in the y-dimension of the environment
%   nbins_env = max number of dispersal bins the env can support (must be >= 2)
%   b = number of offspring per site in the environment
%
% SAVED VARIABLES:
%   E = 2D matrix of patches, marked as viable (1) and nonviable (0)
%   sx = number of sites in the x-dimension of the environment
%   sy = number of sites in the y-dimension of the environment
%   S = number of actual viable sites in the environment
%   xcoord = list of x-coordinates for all patches in E
%   ycoord = list of y-coordinates for all patches in E
%   via_ID = list of viable patches in E, marked by their index number
%   dmap = stucture, mapping from all natal patches to all possible
%      post-dispersal patches
%      dmap{i} is a 2D matrix that maps from row to column patches
%   tmap = 2D matrix the size of E that maps each post-dispersal patch back
%      to a settlement patch (by its element within E) or to death (0)
%   nbins_env = max number of dispersal bins the env can support (must be >= 2)
%   bvec = vector of offspring capacity of each viable site, 1...S

    smax = max(sx,sy); %max environmental distance

% create the environment
    E_via = ones(sy,sx);                 % viable patches in the environment
    E = zeros(sy+2*smax,sx+2*smax);      % full environment
    E(smax+1:smax+sy,smax+1:smax+sx) = E_via;
    via_ID = find(E);                  % ID of viable patches within E
    S = sum(sum(E));                   % total number of viable patches

% get dimensions and coordinates off the environment
    xdim = size(E,2);     % number of patches in x dimension
    ydim = size(E,1);     % number of patches in y dimension
    xcoord = repmat(1:xdim,ydim,1);
    xcoord = xcoord(:);                % vector of x-coordinates for each patch
    ycoord = repmat([1:ydim]',1,xdim);
    ycoord = ycoord(:);                % vector of y-coordinates for each patch

% create map from origin patches to possible post-dispersal patches
    % dmap{1} lists self
    % dmap{2} lists possible dispersal patches at distance 1
    % dmap{j} is a matrix where the rows correspond to each patch in E_via,
    %   and the columns list the ID number in E for all the
    %   patches distance j-1 away from that origin
    dmap{1} = via_ID; % patches in E that are 1, i.e. correspond to E_via
    
    for ii = 1:nbins_env-1
        % for each dispersal distance up to the max do the following:
        for i = 1:S
            % for each viable patch 1 to S, find its corresponding ID within E
            % using via_ID, then use this to get its x and y coordinate within
            % E, and finally use this to get the other patches at distance
            % ~ii away within E based on x and y coordinates.
            tmp(i,:) = find( (xcoord-xcoord(via_ID(i))).^2 + (ycoord-ycoord(via_ID(i))).^2<=ii^2 & (xcoord-xcoord(via_ID(i))).^2 + (ycoord-ycoord(via_ID(i))).^2>(ii-1)^2);            
        end
        dmap{ii+1} = tmp;
        clear tmp
    end
    
% create map from post-dispersal patches to settlement (viable) patches
    tmap = zeros(sy+2*smax,sx+2*smax);
    tmap(smax-sy+1:smax-sy+sy*3,smax-sx+1:smax-sx+sx*3) = repmat(reshape(1:S,sy,sx),3,3);
    % ^ this is a bit hacky for the moment -- it basically replicates the
    % origin patches 3-by-3 and leaves zeros around the edge so that tmap
    % is the same size as E (which we need to be true). This means that
    % both x and y boundaries are wrapped.

% create vector of the number of offspring produced at each site
bvec = b*ones(S,1);


figure(1);
imagesc(E_via);
colormap(flipud(gray))
grid on
saveas(1,strcat(['env_unbounded_sx=' num2str(sx) '_sy=' num2str(sy) '.jpg']))

clear E_via i ii smax tmp xdim ydim 
    
save(strcat(['env_unbounded_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins_env) '.mat']))

