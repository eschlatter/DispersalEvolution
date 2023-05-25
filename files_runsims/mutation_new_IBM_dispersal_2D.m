function [] = mutation_new_IBM_dispersal_2D(gflag,eflag,sx,sy,nbins_env,nbins,nmax,G,b,p,del,pop_init,saveto_filepath)
% written by Allison K. Shaw (ashaw@umn.edu)
% updated May 2022 by E Schlatter
%
% IBM version of Hamilton & May model with spatially-explicit patches in
% 2-D
% Two steps to the dispersal process: individuals go from natal patch to
%   landing patch, then move from landing patch to final patch
%  > S sites in the environment, each supporting K individuals
%  > each adult produces b offspring per generation
%  > each individual is defined by set of dispersal strategies (probability
%      of traveling distance x)
%  > each dispersing individual has a probability p of surviving dispersal
%      (regardless of distance traveled)
%  > offspring inherit dispersal strategy from parent with small mutation
%  > individuals whose landing patch is not viable or who disperse beyond
%     the edge of the world die...except if they navigate back to a close
%     viable patches
%
% INPUTS:
%   gflag = whether (1) or not (0) to show plots during simulation
%   eflag = which environment to use: 1=unbounded, 2=bounded, 3=reef,
%   4=heterogeneous
%   sx = number of sites in the x-dimension of the environment
%   sy = number of sites in the y-dimension of the environment
%   nbins_env = max number of dispersal bins the env can support (must be >= 2)
%   nbins = how many dispersal bins to actually use
%   nmax = maximum larval navigation distance (behavior)
%   G = number of total generations to simulate
%   b = offspring produced per individual
%   p = probability of surviving dispersal
%   del = fraction of dispersal probability to move during mutation
% 

rng('shuffle') % seed the random number generator from computer clock

K = 1;         % carrying capacity per patch

%saveto_filepath = '../output_simulations/20221101';

if nbins < 2; error('nbins_use must be at least 2'); end
if nbins > nbins_env; error('nbins_env must be bigger than nbins'); end

%-----LOAD-ENVIRONMENT----------------------------------------------------%
    if eflag==1
        load(strcat(['../output_environments/env_unbounded_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins_env) '.mat']))
    elseif eflag==2
        load(strcat(['../test_output_environments/env_bounded_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins_env) '_nmax=' num2str(nmax) '.mat']))
    elseif eflag==3
        load(strcat(['../output_environments/env_reef_nbins=' num2str(nbins_env) '_nmax=' num2str(nmax) 'dmap_bmin=' num2str(b(1)) '_bmax=' num2str(b(2)) '.mat']))
    elseif eflag==4
        load(strcat(['../output_environments/env_bounded_het_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins_env) '_nmax=' num2str(nmax) '_bmin=' num2str(b(1)) '_bmax=' num2str(b(2)) '.mat']))
    end
%-----LOAD-ENVIRONMENT----------------------------------------------------%


%-----INITALIZATION-------------------------------------------------------%
    N0 = round(K*S); % initial number of individuals - at capacity

    %check if a specified initial population matrix matches the rest of the
    %parameters
     % if pop_init is the right dimensions, 
     % REMOVED THIS % and all the displacement strategies add up to 1  && isequal(sum(pop_init(:,1:nbins),2),ones(1,N0)') 
      % and contains only viable patches
    if isequal(size(pop_init),[N0,nbins+1]) && sum(ismember(pop_init(:,nbins+1),via_ID))==N0 
        pop=pop_init;
    else
        pop_init=0;
        sprintf('pop_init incorrectly specified; using default starting population')
    end

    if pop_init==0 % if no initial population matrix specified
        % create matrix to hold information for all individuals
        pop = zeros(N0,nbins+1);
        pop(:,1:nbins) = (1/nbins)*ones(N0,nbins); % assign starting dispersal bin values (set equal prob of each distance)
        pop(:,nbins+1) = repmat(via_ID,1,K); % assign evenly to starting patches (patch named by absolute index (in xcoord,etc)
    end

    % matrices to record population parameters in, over time
    dtime_avg = zeros(G,nbins); % population mean of dispersal parameters
    dtime_std = zeros(G,nbins); % pop. standard deviation of dispersal parammeters

    mortality_cost = zeros(G,nbins);
    kincomp_cost = zeros(G,nbins);

%     % matrix to record survival rates over time
%     fitness = zeros(G,3);

    % matrices to hold displacement, dispersal, and recruitment kernels
    % from individual dispersal events
%     nbins_plus = nbins+nmax;
%     kernel_displacement = zeros(G,nbins_plus);
%     kernel_dispersal = zeros(G,nbins_plus);
%     kernel_recruitment = zeros(G,nbins_plus);

%     % matrices to hold parent-level larval survival rates
%     parent_kernel_mean = zeros(length(via_ID),G); %mean of the parental displacement strategy
%     parent_surv_dispersal = zeros(length(via_ID),G); %larvae that make it to a habitable site (by displacement or navigation)
%     parent_surv_recruitment = zeros(length(via_ID),G); %larvae that actually recruit

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
%     if(size(pop,1)==length(via_ID))
%         parent_kernel_mean(:,g) = pop(:,1:nbins)*(0:(nbins-1))'; %mean of each parent's displacement strategy
%     end

    %-----REPRODUCTION-----%
        % each individual produces b offspring with same patch and dispersal
        % parameters as parent
        off = [];
        for j = 1:size(pop,1)
            off = [off; repmat(pop(j,:),[patches(pop(j,nbins+1)),1])];
        end

        off(1,nbins+2:nbins+5)=[0 0 0 0]; %add columns to off to store more info
        Noff = size(off,1); % total number of offspring
    
        % off: columns
        %  1:nbins > displacement kernel
        %  nbins+1 > starting patch
        %  nbins+2 > displacement patch (or 0, if died during displacement)
        %  nbins+3 > dispersal patch (or 0, if died during navigation)
        %  nbins+4 > recruitment patch (or 0, if died during competition)
        %  nbins+5 > displacement distance

%         fitness(g,1) = Noff; % store total number of offspring produced
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

%     % store displacement distances
%     kernel_displacement(g,:) = sum(d_ind == 1:nbins_plus);

    clear ind val j

    off(:,nbins+5) = d_ind; % store displacement distances in off
    % d_ind starts at 1 for offspring who stay home

    ind = find(d_ind==1); % which offspring didn't leave natal patch during displacement
    off(ind,nbins+2) = off(ind,nbins+1); % their patch post-displacement is same as patch pre-displacement
    clear ind

    for j = 2:nbins     % for each of the possible dispersal distances
        ind = find(d_ind==j); % find all offspring that drew j-1 dispersal distance
        for i = ind'
            if surv(i)==1 % for offspring that survived dispersal
                x = find(dists(:,off(i,nbins+1))==(j-1)); %possible patches to disperse to
                y = randi(length(x)); % pick one at random
                off(i,nbins+2)=x(y); % save landing patch
                clear x y
            end
        end

    end

    % remove offspring who died during dispersal here
%     died = off(:,nbins+2)==0; % any offspring not assigned a patch above died during displacement
%     off(died,:) = []; % remove offspring who died during displacement
% %     fitness(g,2) = size(off,1); % record number of offspring left after displacement mortality
%     clear died

    % navigation
    for i = 1:size(off,1) %for each remaining larva
        if surv(i)==1 %only larvae that didn't die during displacement
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
    end

%     % find offspring who didn't settle anywhere and remove
%     died = off(:,nbins+3)==0;
%     off(died,:) = []; % remove offspring who died during navigation
% %     fitness(g,3) = size(off,1); % record number of offspring left after navigation
%     clear died

%     % store number of offspring who survived dispersal from each original patch
%     for patch = 1:length(via_ID)
%         parent_surv_dispersal(patch,g)=sum(off(:,nbins+1)==via_ID(patch));
%     end
%     clear patch

%     % calculate and store dispersal distances
%     idx = sub2ind(size(dists),off(:,nbins+1),off(:,nbins+3));
%     dispersal_distances = dists(idx);
%     clear idx
%     dispersal_distances = ((xcoord(off(:,nbins+1))-xcoord(off(:,nbins+3))).^2 + (ycoord(off(:,nbins+1))-ycoord(off(:,nbins+3))).^2).^0.5;
%     dispersal_distances = floor(dispersal_distances);
%     kernel_dispersal(g,:) = sum(dispersal_distances == 0:(nbins_plus-1));

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
    %Noffs = sum(off(:,nbins+3)==via_ID'); % number of offspring per patch
    %fullind = find(Noffs>K);  % index of overcrowded patches
    for i = 1:length(via_ID) % loop over each patch
        patch_ind = find(off(:,nbins+3)==via_ID(i)); % index of offspring in patch
        if(~isempty(patch_ind))
            n_settlers = min(K,size(patch_ind,1));
            off(patch_ind(1:n_settlers),nbins+4) = off(patch_ind(1:n_settlers),nbins+3); %choose K to survive (gets a patch in column nbins+4 = recruitment column)
        end
    end

    clear fullind patch_ind i
    % save recruited offspring as new population
    remaining = find(off(:,nbins+4)~=0);
    pop = off(remaining,[1:nbins,nbins+3]);
    %-----COMPETITION-----%

    % store number of offspring who survived recruitment from each original patch
%     for patch = 1:length(via_ID)
%         parent_surv_recruitment(patch,g)=sum(off(:,nbins+1)==via_ID(patch));
%     end
%     clear patch

%     % calculate and store recruitment distances
%     idx = sub2ind(size(dists),off(:,nbins+1),off(:,nbins+3));
%     recruitment_distances = dists(idx);
%     clear idx
%     recruitment_distances = ((xcoord(off(:,nbins+1))-xcoord(off(:,nbins+3))).^2 + (ycoord(off(:,nbins+1))-ycoord(off(:,nbins+3))).^2).^0.5;
%     recruitment_distances = floor(recruitment_distances);
%     kernel_recruitment(g,:) = sum(recruitment_distances == 0:(nbins_plus-1));

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

    %-----CALCULATE OUTPUTS-----%

    % direct cost (proportion of larvae that die because they can't reach
    % suitable habitat)
    for i = 1:nbins
        off_dist = off(off(:,nbins+5) == i,:);
        died_in_nav = sum(off_dist(:,nbins+3)==0);
        mortality_cost(g,i) = died_in_nav/size(off_dist,1);
    end
    clear off_dist died_in_nav

    
    % indirect cost (expected number of siblings encountered by a
    % potentially-recruiting larva at its destination site, divided by total
    % number of competitors at that site)
    for i = 1:nbins
        % choose the offspring that displaced the focal distance AND survived to compete
        off_dist = off(off(:,nbins+5) == i,:);
        off_dist = off_dist(off_dist(:,nbins+3) ~= 0,:);

        kincost_i=zeros(size(off_dist,1),1);
        % pick each individual in off_dist. Find its destination site. Then
        % count the total number of individuals (from off) in that
        % destination site, and the number that share the focal
        % individual's origin site.
        for j = 1:size(off_dist,1)
            origin = off_dist(j,nbins+1);
            destination = off_dist(j,nbins+3);
            %n_competitors = sum(off(:,nbins+3)==destination);
            sibs = sum((off(:,nbins+1)==origin).*(off(:,nbins+3)==destination))-1;
            %kincost_i = [kincost_i sibs/n_competitors];
            kincost_i(j) = sibs/b;
        end
        clear origin destination
        
        kincomp_cost(g,i) = mean(kincost_i);
        clear off_dist
    end

%-----CALCULATE OUTPUTS-----%

end % generation loop

toc

    % save output as mat file
    if eflag==1
        save(strcat([saveto_filepath '/IBM_unbounded_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins) '_nmax=' num2str(nmax) '_del=' num2str(del) '_b=' num2str(b) '_p=' num2str(p) '.mat']))
    elseif eflag==2
        save(strcat([saveto_filepath '/IBM_bounded_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins) '_nmax=' num2str(nmax) '_del=' num2str(del) '_b=' num2str(b) '_p=' num2str(p) '.mat']))
    elseif eflag==3
        save(strcat([saveto_filepath '/IBM_reef_nbins=' num2str(nbins) '_nmax=' num2str(nmax) '_bmin=' num2str(bmin) '_bmax=' num2str(bmax) '_del=' num2str(del) '_p=' num2str(p) '.mat']))
    elseif eflag==4
        save(strcat([saveto_filepath '/IBM_bounded_het_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins_env) '_nmax=' num2str(nmax) '_bmin=' num2str(bmin) '_bmax=' num2str(bmax) '_del=' num2str(del) '_p=' num2str(p) '.mat']))
    end

%-----SIMULATE------------------------------------------------------------%