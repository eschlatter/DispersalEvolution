function [] = create_env_bounded(sx,sy,nbins_env,nmax,b)
% written by Allison K. Shaw (ashaw@umn.edu)
% updated: June 2018
%
% USAGE: [] = create_env_bounded(sx,sy,nbins_env,nmax,b)
%
% Create the environment to use as a backdrop for the dispersal IBM.
% Here, the environment is a set of viable patches, surrounded by
%   nonviable patches.
% Only need to run this once - it creates a mat file that can be uploaded
%   later within IBM code.
%
% INPUTS:
%   sx = number of sites in the x-dimension of the environment
%   sy = number of sites in the y-dimension of the environment
%   nbins_env = max number of dispersal bins the env can support (must be >= 2)
%   nmax = maximum larval navigation distance (behavior)
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
%   nmax = maximum larval navigation distance (behavior)
%   bvec = vector of offspring capacity of each viable site, 1...S
% %% dists = array of distances between each site

    stmp = nbins_env+1; % make sure environment can hold all max dispersal distances

% create the environment
    E_via = ones(sy,sx);                 % viable patches in the environment
    E = zeros(sy+2*stmp,sx+2*stmp);      % full environment
    E(stmp+1:stmp+sy,stmp+1:stmp+sx) = E_via;
    via_ID = find(E);                  % ID of viable patches within E
    S = sum(sum(E));                   % total number of viable patches

% get dimensions and coordinates off the environment
    xdim = size(E,2);     % number of patches in x dimension
    ydim = size(E,1);     % number of patches in y dimension
    xcoord = repmat(1:xdim,ydim,1);
    xcoord = xcoord(:);                % vector of x-coordinates for each patch
    ycoord = repmat([1:ydim]',1,xdim);
    ycoord = ycoord(:);                % vector of y-coordinates for each patch

% %% create array of distances
    dists=zeros(length(xcoord),length(xcoord));
    for ii = 1:length(xcoord)
        for i = 1:length(xcoord)
            dists(ii,i)=((xcoord(ii)-xcoord(i))^2+(ycoord(ii)-ycoord(i))^2)^.5;
        end
    end
    dists=floor(dists);
    
% % create map from origin patches to possible post-dispersal patches
%     % dmap{1} lists self
%     % dmap{2} lists possible dispersal patches at distance 1
%     % dmap{j} is a matrix where the rows correspond to each patch in E_via,
%     %   and the columns list the ID number in E for all the
%     %   patches distance j-1 away from that origin
%     dmap{1} = via_ID; % patches in E that are 1, i.e. correspond to E_via
%     
%     for ii = 1:nbins_env-1
%         % for each dispersal distance up to the max do the following:
%         for i = 1:S
%             % for each viable patch 1 to S, find its corresponding ID within E
%             % using via_ID, then use this to get its x and y coordinate within
%             % E, and finally use this to get the other patches at distance
%             % ~ii away within E based on x and y coordinates.
%             tmp(i,:) = find( (xcoord-xcoord(via_ID(i))).^2 + (ycoord-ycoord(via_ID(i))).^2<=ii^2 & (xcoord-xcoord(via_ID(i))).^2 + (ycoord-ycoord(via_ID(i))).^2>(ii-1)^2);            
%         end
%         dmap{ii+1} = tmp;
%         clear tmp
%     end
%   
% % create map from post-dispersal patches to settlement (viable) patches
%     tmap = zeros(sy+2*stmp,sx+2*stmp);
%     % start by mapping viable patches back to themselves
%     tmap(find(E))=1:S;
%     %-----LARVAL-BEHAVIOR-----%   
%         % x and y coordinate of viable patches
%         xcoord_via = xcoord(tmap>0);
%         ycoord_via = ycoord(tmap>0);
% 
%         for ii = 1:nmax % loop up to maximum larval navigation distance
% 
%             % index of non-viable patches
%             ind = find(tmap==0);
% 
%             for i = 1:length(ind)  %for each empty patch
% 
%                 % find viable patches, if any within distance ii
%                 tmp = find( (xcoord_via-xcoord(ind(i))).^2 + (ycoord_via-ycoord(ind(i))).^2<=ii^2);
% 
%                 if tmp % if there are any
%                     y = randi(length(tmp)); % pick one at random
%                     tmap(ind(i)) = tmp(y); % save it
%                 end
% 
%             end
%         end
%     %-----LARVAL-BEHAVIOR-----%   

% create vector of the number of offspring produced at each site
bvec = b*ones(S,1);

% %% create vector patches, with birth rate of each site (viable), or 0
% (nonviable)
patches = zeros(length(xcoord),1);
patches(via_ID)=bvec; %added to create_env_bounded

clear E_via i ii ind stmp tmp xcoord_via ycoord_via xdim ydim y

save(strcat(['../test_output_environments/env_bounded_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins_env) '_nmax=' num2str(nmax) '.mat']),'-v7.3')

figure(1);
imagesc(E);
colormap(flipud(gray))
grid on
saveas(1,strcat(['../test_output_environments/env_bounded_sx=' num2str(sx) '_sy=' num2str(sy) '.jpg']))
