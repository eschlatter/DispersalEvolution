%% Parameters
clear
%BIOLOGICAL PARAMETERS
gflag = 0;     % whether (1) or not (0) to show plots during simulation
b = 10;        % offspring produced per individual
del = 0.001;   % fraction of dispersal probability to move during mutation
p = 1;      % probability of surviving dispersal

%ENVIRONMENT PARAMETERS
nbins_env = 5; % max number of dispersal bins used when creating env (>=2)
nbins = 5; % how many dispersal bins to actually use (2<=nbins_use<=nbins_max)
eflag = 2; % which environment to use: 1=unbounded, 2=bounded, 3=reef
S = 32^2;      % number of sites in the environment
sx = 2; % number of sites in the x-dimension of the environment
sy = S/sx;
    
G = 5; % number of total generations to simulate

% larval navigation distance of 0, 1, 2 or 3
nmax = 1;  % maximum larval recruitment distance (behavior)

rng('shuffle') % seed the random number generator from computer clock

K = 1;         % carrying capacity per patch

saveto_filepath = '../output_simulations/20220516_test';

if nbins < 2; error('nbins_use must be at least 2'); end
if nbins > nbins_env; error('nbins_env must be bigger than nbins'); end

%-----LOAD-ENVIRONMENT----------------------------------------------------%
    if eflag==1
        load(strcat(['../output_environments/env_unbounded_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins_env) '.mat']))
    elseif eflag==2
        load('C:\Users\eschlatter\Dropbox\DispersalEvolution\test_output_environments\env_bounded_sx=2_sy=512_nbins=5_nmax=1.mat')
        %load(strcat(['../test_output_environments/env_bounded_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins_env) '_nmax=' num2str(nmax) '.mat']))
    elseif eflag==3
        load(strcat(['../output_environments/env_reef_nbins=' num2str(nbins_env) '_nmax=' num2str(nmax) 'dmap_bmin=' num2str(b(1)) '_bmax=' num2str(b(2)) '.mat']))
    elseif eflag==4
        load(strcat(['../output_environments/env_bounded_het_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins_env) '_nmax=' num2str(nmax) '_bmin=' num2str(b(1)) '_bmax=' num2str(b(2)) '.mat']))
    end
%-----LOAD-ENVIRONMENT----------------------------------------------------%

%% Simulation

%-----INITALIZATION-------------------------------------------------------%
    N0 = round(K*S); % initial number of individuals - at capacity

    % create matrix to hold information for all individuals
    pop = zeros(N0,nbins+1);
    pop(:,1:(nmax+2)) = (1/(nmax+2))*ones(N0,nmax+2); % assign starting dispersal bin values (set equal prob of each distance up to navigation distance + 1)
    pop(:,nbins+1) = repmat(via_ID,1,K); % assign evenly to starting patches (patch named by absolute index (in xcoord,etc)
    
    % matrices to record population parameters in, over time
    dtime_avg = zeros(G,nbins); % population mean of dispersal parameters
    dtime_std = zeros(G,nbins); % pop. standard deviation of dispersal parammeters

    % matrix to record survival rates over time
    fitness = zeros(G,3);

    % matrices to hold displacement, dispersal, and recruitment kernels
    % from individual dispersal events
    nbins_plus = nbins+nmax;
    kernel_displacement = zeros(G,nbins_plus);
    kernel_dispersal = zeros(G,nbins_plus);
    kernel_recruitment = zeros(G,nbins_plus);

%%%%%%%%%%%%%    
    % matrices to hold parent-level larval survival rates
    parent_surv_dispersal = zeros(length(via_ID),G);
    parent_surv_recruitment = zeros(length(via_ID),G);

    % if set to display graphics, display first figure and wait for keystrike
    if gflag==1
        figure(1); clf
        plot(pop(:,1:nbins)','o')
        ylim([0 1])
        xlabel('dispersal bin')
        ylabel('probability')
        pause(0.1)
    end
%-----INITALIZATION-------------------------------------------------------%

tic

%-----SIMULATE------------------------------------------------------------%
g=0; % this counts the number of generations that have passed

while g<G && size(pop,1)>0 % loop over generations (only while population not extinct)
    g = g+1 % step generation forward

    % record statistics
    dtime_avg(g,:) = mean(pop(:,1:nbins),1);   % average dispersal params
    dtime_std(g,:) = std(pop(:,1:nbins),[],1); % std. of dispersal params

    %-----REPRODUCTION-----%
        % each individual produces b offspring with same patch and dispersal
        % parameters as parent
        off = [];
        for j = 1:size(pop,1)
            off = [off; repmat(pop(j,:),[patches(pop(j,nbins+1)),1])];
        end
        off(1,nbins+2:nbins+3)=[0 0]; %add two more columns to off (will hold displacement and dispersal patches)
        Noff = size(off,1); % total number of offspring
        fitness(g,1) = Noff; % store total number of offspring produced
    %-----REPRODUCTION-----%


    %-----MUTATION-AND-DISPERSAL-----%
    srand = rand(Noff,1); % random numbers to use for survival probabilities
    surv = srand<p;        % which offspring survive (1) or not (0), if they attempt to disperse
    drand = rand(Noff,1); % random numbers to use for dispersal probabilities
    d_ind = zeros(Noff,1); % to hold the index of dispersal distance traveled (1=stay, 2=distance 1, 3=distance 2, etc)

    for j = 1:Noff

        % mutate dispersal strategy slightly by taking or adding del from dispersal bin
        % if del amount isn't remaining in dispersal bin, take everything
        % (i.e. bound dispersal probability at zero)
        ind = randperm(nbins,2);         % select two bins at random
        val = min(off(j,ind(1)),del);          % take del of prob, if available
        off(j,ind) = off(j,ind) + [-val +val]; % move del from first to second

        % find dispersal bin corresponding to random dispersal distance
        % cumsum is the cumulative sum across the disperal bins -- i.e. the
        %   cumulative probability distribution
        % this finds the first bin where cumsum exceeds the random
        %   dispersal number generated
        % e.g. if off(j,1:nbins) = [0.1 0.8 0.1] this is 10%, 80%, and
        %      10% probability of traveling distances 0, 1, and 2
        %      and cumsum(off(j,1:nbins)) = [0.1 0.9 1]
        %      if drand(j) = 0.2 then d_ind(j) would be 2, that is the
        %      individual travels distance 1 during dispersal
        d_ind(j) = find(cumsum(off(j,1:nbins),2)>drand(j),1,'first');

    end

    % store displacement distances
    kernel_displacement(g,:) = sum(d_ind == 1:nbins_plus);

    clear ind val j

    ind = find(d_ind==1); % which offspring didn't leave natal patch during displacement
    off(ind,nbins+2) = off(ind,nbins+1); % their patch post-displacement is same as patch pre-displacement
    clear ind

    for j = 2:nbins     % for each of the possible dispersal distances
        ind = find(d_ind==j); % find all offspring that drew j-1 dipersal distance
        for i = ind'
            % for offspring that survived dispersal
            if surv(i)==1
                x = find(dists(:,off(i,nbins+1))==(j-1)); %possible patches to disperse to
                y = randi(length(x)); % pick one at random
                off(i,nbins+2)=x(y); % save landing patch
% %%                Note that we've switched to Euclidean distance here (vs
%                       the von Neumann distance that comes from dmap)
            end
        end

    end

    % remove offspring who died during dispersal here
    died = off(:,nbins+2)==0; % any offspring not assigned a patch above died during displacement
    off(died,:) = []; % remove offspring who died during displacement
    fitness(g,2) = size(off,1); % record number of offspring left after displacement mortality
    clear died

    % navigation
    for i = 1:size(off,1) %for each remaining larva
        if patches(off(i,nbins+2))==0 %if the larva has displaced to an uninhabitable patch
            x = find(dists(:,off(i,nbins+2))<=nmax); %find patches within navigation distance
            x_hab = x(patches(x)~=0); %restrict to habitable patches
            if ~isempty(x_hab) %if there are any habitable patches
                y = randi(length(x_hab)); %pick one at random
                off(i,nbins+3) = x_hab(y); %save settlement patch
            end
        else % if the larva has displaced to a habitable patch
            off(i,nbins+3) = off(i,nbins+2); %stay there
        end
    end

    % find offspring who didn't settle anywhere and remove
    died = off(:,nbins+3)==0;
    off(died,:) = []; % remove offspring who died during navigation
    fitness(g,3) = size(off,1); % record number of offspring left after navigation
    clear died

%%%%%%%%%%%%%%%%%    
    % store number of offspring who survived dispersal from each original patch
    for patch = 1:length(via_ID)
        parent_surv_dispersal(patch,g)=sum(off(:,6)==via_ID(patch));
    end

    % calculate and store dispersal distances
    idx = sub2ind(size(dists),off(:,nbins+1),off(:,nbins+3));
    dispersal_distances = dists(idx);
    clear idx
% %%     dispersal_distances = ((xcoord(off(:,nbins+1))-xcoord(off(:,nbins+3))).^2 + (ycoord(off(:,nbins+1))-ycoord(off(:,nbins+3))).^2).^0.5;
% %%     dispersal_distances = floor(dispersal_distances);
    kernel_dispersal(g,:) = sum(dispersal_distances == 0:(nbins_plus-1));

    clear srand surv drand d_ind ind i j
    %-----MUTATION-AND-DISPERSAL-----%


    %-----COMPETITION-----%
    % randomly reorder all offspring (to avoid spurious patterns in next
    % step)
    Noff = size(off,1); % total number of offspring
    xind = randperm(Noff);
    off = off(xind,:);

    % only allow K offspring per patch to survive
    % %%(checked - all remaining offspring are in viable patches)
    Noffs = sum(off(:,nbins+3)==via_ID'); % number of offspring per patch
    fullind = find(Noffs>K);  % index of overcrowded patches
    for i = 1:length(fullind) % loop over each overcrowded patch
        patch_ind = find(off(:,nbins+3)==via_ID(fullind(i))); % index of offspring in patch
% %%    patch_ind = find(off(:,nbins+1)==via_ID(fullind(i))); % !! typo
        off(patch_ind(K+1:end),:) = []; % kill all but K of these
    end
    clear fullind patch_ind i
    % kill off all adults and just save offspring as new population
    pop = off(:,[1:nbins,nbins+3]);
    %-----COMPETITION-----%

%%%%%%%%%%%%%%%%%    
    % store number of offspring who survived recruitment from each original patch
    for patch = 1:length(via_ID)
        parent_surv_recruitment(patch,g)=sum(off(:,6)==via_ID(patch));
    end

    % calculate and store recruitment distances
    idx = sub2ind(size(dists),off(:,nbins+1),off(:,nbins+3));
    recruitment_distances = dists(idx);
    clear idx
% %%    recruitment_distances = ((xcoord(off(:,nbins+1))-xcoord(off(:,nbins+3))).^2 + (ycoord(off(:,nbins+1))-ycoord(off(:,nbins+3))).^2).^0.5;
% %%    recruitment_distances = floor(recruitment_distances);
    kernel_recruitment(g,:) = sum(recruitment_distances == 0:(nbins_plus-1));

    % if set to display graphics, update figure 1
    if gflag==1
        N = size(pop,1);  % count number of individuals present
        figure(1); clf
        plot(pop(:,1:nbins)','o')
        ylim([0 1])
        xlabel('dispersal bin')
        ylabel('probability')
        title(strcat(['generation = ' num2str(g)]))
        pause(0.1)
    end

end % generation loop

%-----SIMULATE------------------------------------------------------------%




%%%%%%%%%%%%%
% plot # of offspring that successfully disperse
figure
clf
mean_parent_surv_dispersal = mean(parent_surv_dispersal);
sd_parent_surv_dispersal = std(parent_surv_dispersal);
errorbar(mean_parent_surv_dispersal,sd_parent_surv_dispersal)
title(sprintf('Per-parent dispersal survival rate, nmax = %g',nmax));
xlabel('Generation')
ylabel('Successfully dispersing offspring per parent (mean+/-sd)')

% plot # of offspring that successfully recruit
figure
clf
mean_parent_surv_recruitment = mean(parent_surv_recruitment);
sd_parent_surv_recruitment = std(parent_surv_recruitment);
errorbar(mean_parent_surv_recruitment,sd_parent_surv_recruitment)
title(sprintf('Per-parent recruitment rate, nmax = %g',nmax));
xlabel('Generation')
ylabel('Successfully recruiting offspring per parent (mean+/-sd)')