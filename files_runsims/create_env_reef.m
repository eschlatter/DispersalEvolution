function [] = create_env_reef(nbins_env,nmax,bmin,bmax)
% written by Allison K. Shaw (ashaw@umn.edu)
% updated: March 2016
%
% USAGE: [] = create_env_reef(nbins_env,nmax)
%
% Create the environment to use as a backdrop for the dispersal IBM.
% Here, the environment is a copy of the Belize reef with viable patches,
%   surrounded by nonviable patches.
% Only need to run this once - it creates a mat file that can be uploaded
%   later within IBM code.
%
% INPUTS:
%   nbins_env = max number of dispersal bins the env can support (must be >= 2)
%   nmax = maximum larval navigation distance (behavior)
%   bmin = minimum offspring per site in the environment
%   bmax = maximum offspring per site in the environment
%
% SAVED VARIABLES:
%   E = 2D matrix of patches, marked as viable (1) and nonviable (0)
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
%   nmax = maximum larval navigation distance (behavior)



% Import Belize data from ascii file
%    Attached is an ACII file of 0s and 1s that is based on the raster file
%    of 1 km square cells. Unsuitable habitat is 0 and forereef is 1. The
%    cell size says 1,000 because the units are meters.
    % import data using space as delimiter and skipping 6 header rows
    A = importdata('trimmed_ascii.txt',' ',6);
    E_tmp = A.data; % 626 rows by 874 columns
    % replace -9999 with 0s (just first row for some reason)
    i = find(E_tmp==-9999);
    E_tmp(i) = 0;
    % trim data down
    E_tmp = E_tmp(392:626,700:797);
    % add zero padding
    E_tmp = [zeros(1,size(E_tmp,2)); E_tmp; zeros(1,size(E_tmp,2))];

% calculate x and y dimensions of the viable patches
    [sy,sx] = size(E_tmp);
    smax = max(sx,sy);

% create the environment
    E = zeros(sy+2*smax,sx+2*smax);      % full environment
    E(smax+1:smax+sy,smax+1:smax+sx) = E_tmp;
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
    % start by mapping viable patches back to themselves
    tmap(find(E))=1:S;
    %-----LARVAL-BEHAVIOR-----%   
        % x and y coordinate of viable patches
        xcoord_via = xcoord(tmap>0);
        ycoord_via = ycoord(tmap>0);

        for ii = 1:nmax % loop up to maximum larval navigation distance

            % index of non-viable patches
            ind = find(tmap==0);

            for i = 1:length(ind)  %for each empty patch

                % find viable patches, if any within distance ii
                tmp = find( (xcoord_via-xcoord(ind(i))).^2 + (ycoord_via-ycoord(ind(i))).^2<=ii^2);

                if tmp % if there are any
                    y = randi(length(tmp)); % pick one at random
                    tmap(ind(i)) = tmp(y); % save it
                end

            end
        end
    %-----LARVAL-BEHAVIOR-----%   

% create vector of the number of offspring produced at each site
% (a random integer between bmin and bmax)
bvec = randi([bmin bmax],S,1);

clear A E_tmp i ii ind smax sx sy tmp xcoord_via ycoord_via xdim ydim y

save(strcat(['env_reef_nbins=' num2str(nbins_env) '_nmax=' num2str(nmax) '_bmin=' num2str(bmin) '_bmax=' num2str(bmax) '.mat']))

figure(1);
Etmp = E;
Etmp(find(Etmp))=bvec;
imagesc(Etmp);
colormap(flipud(gray))
grid on
saveas(1,strcat(['env_reef_bmin=' num2str(bmin) '_bmax=' num2str(bmax) '.jpg']))
clear Etmp
