clear all; clc; close all
% written by Allison K. Shaw (ashaw@umn.edu)
% updated: October 2018
% get counts for Fig 5

% PARAMETERS
nbins_max = 31; % how many dispersal bins to actually use (2<=nbins_use<=nbins_max)
nbins = 31; % how many dispersal bins to actually use (2<=nbins_use<=nbins_max)
nmax = 0; % maximum larval navigation distance (behavior)
b0 = 10;        % offspring produced per individual
bmin = b0; % minimum offspring per site in the environment
bmax = b0; % maximum offspring per site in the environment
del = 0.001;   % fraction of dispersal probability to move during mutation
p = 1; % probability of surviving dispersal

counts_1 = NaN(nbins,3);
counts_2 = NaN(nbins,3);
counts_4 = NaN(nbins,3);

rng('shuffle') % seed the random number generator from computer clock
for K = [1 2 4]
    
    if K==1
        load(strcat(['../output_simulations/IBM_reef_nbins=' num2str(nbins) '_nmax=' num2str(nmax) '_bmin=' num2str(bmin) '_bmax=' num2str(bmax) '_del=' num2str(del) '_p=' num2str(p) '.mat']))
    else
        load(strcat(['../output_simulations/IBM_reef_nbins=' num2str(nbins) '_nmax=' num2str(nmax) '_bmin=' num2str(bmin) '_bmax=' num2str(bmax) '_del=' num2str(del) '_p=' num2str(p) '_K=' num2str(K) '.mat']))
    end
    
    clear E G N N0 Noff Noffs b dtime_avg dtime_std eflag
    clear g gflag nbins_env pat_sett
    clear off sx sy x xind y
    
    %%-------------FROM IBM CODE VERBATIM
            %-----REPRODUCTION-----%
            % each individual produces b offspring with same patch and dispersal
            % parameters as parent
            off = [];
            for j = 1:size(pop,1)
                off = [off; repmat(pop(j,:),[bvec(pop(j,nbins+1)),1])];
            end
            Noff = size(off,1); % total number of offspring
            %-----REPRODUCTION-----%


            %-----MUTATION-AND-DISPERSAL-----%
            srand = rand(Noff,1); % random numbers to use for survival probabilities
            surv = srand<p;        % which offspring survive (1) or not (0), if they attempt to disperse
            drand = rand(Noff,1); % random numbers to use for dispersal probabilities
            d_ind = zeros(Noff,1); % to hold the index of dispersal distance traveled (1=stay, 2=distance 1, 3=distance 2, etc)

            for j = 1:Noff;

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
            clear ind val j
    %%-------------FROM IBM CODE VERBATIM
    
    
% COUNT VERSION 1: PRE DISPERSAL
if     K==1
    counts_1(:,1) = hist(d_ind,1:nbins_max);
elseif K==2
    counts_2(:,1) = hist(d_ind,1:nbins_max);
elseif K==4
    counts_4(:,1) = hist(d_ind,1:nbins_max);
end

    
    %%-------------FROM IBM CODE VERBATIM
           pat_disp = zeros(Noff,1); % to hold the patch where each offspring lands after dispersal

            ind = find(d_ind==1); % which offspring did not disperse
            pat_disp(ind) = dmap{1}(off(ind,nbins+1)); % map natal patch from dmap
            clear ind

            % for each of the possible dispersal distances
            for j = 2:nbins

                ind = find(d_ind==j); % find all offspring that drew j-1 dipersal distance

                for i = ind'
                    % for offspring that survived dispersal
                    if surv(i)==1;
                        x = dmap{j}(off(i,nbins+1),:); % possible patches to disperse to
                        y = randi(length(x)); % pick one at random
                        pat_disp(i) = x(y); % save landing patch
                    end
                end

            end

            % remove offspring who died during dispersal here
            off(pat_disp==0,:) = []; % remove offspring who died during dispersal
            pat_disp(pat_disp==0) = [];

            pat_sett = tmap(pat_disp); % where offspring settle after dispersal
    %%-------------FROM IBM CODE VERBATIM

    
    % remove offspring who died during dispersal here
    d_ind(pat_disp==0) = [];

    
    % record starting and ending patches
    p_start = off(:,nbins_max+1);
    p_end = pat_sett;
    % remove offspring who didn't settle anywhere
    p_start(pat_sett==0) = [];
    p_end(pat_sett==0) = [];

    
        
    %%-------------FROM IBM CODE VERBATIM
           off(:,nbins+1) = pat_sett; % save new locations

            % find offspring who didn't settle anywhere and remove
            off(off(:,nbins+1)==0,:) = []; % remove offspring who died during dispersal

            clear srand surv drand ind i j pat_disp
            %-----MUTATION-AND-DISPERSAL-----%
    %%-------------FROM IBM CODE VERBATIM
    
    
    
    d_ind(pat_sett==0,:) = []; % remove offspring who died during dispersal

    
    % COUNT VERSION 2: POST SURVIVAL, PRE COMPETITION
    if     K==1
        counts_1(:,2) = hist(d_ind,1:nbins_max);
    elseif K==2
        counts_2(:,2) = hist(d_ind,1:nbins_max);
    elseif K==4
        counts_4(:,2) = hist(d_ind,1:nbins_max);
    end
    
    
    %%-------------FROM IBM CODE EDITIED

        %-----COMPETITION-----%
        % randomly reorder all offspring (to avoid spurious patterns in next
        % step)
        Noff = size(off,1); % total number of offspring
        xind = randperm(Noff);
        off = off(xind,:);

        % only allow K offspring per patch to survive
        Noffs = hist(off(:,nbins+1),1:S); % number of offspring per patch
        fullind = find(Noffs>K);  % index of overcrowded patches
        for i = 1:length(fullind) % loop over each overcrowded patch
            patch_ind = find(off(:,nbins+1)==fullind(i)); % index of offspring in patch
            off(patch_ind(K+1:end),:) = []; % kill all but K of these
        d_ind(patch_ind(K+1:end),:) = [];                                  %******************
        end
        clear fullind patch_ind i
        % kill off all adults and just save offspring as new population
        pop = off;
        %-----COMPETITION-----%

    %%-------------FROM IBM CODE EDITIED
    
    
    % convert starting and ending patch to their ID within matrix E
    p_start = via_ID(p_start);
    p_end =  via_ID(p_end);
    p_start_x = xcoord(p_start);
    p_start_y = ycoord(p_start);
    p_end_x = xcoord(p_end);
    p_end_y = ycoord(p_end);

    dist = sqrt((p_start_x - p_end_x).^2 + (p_start_y - p_end_y).^2);
    
    % COUNT VERSION 3: POST SURVIVAL, POST COMPETITION
    if     K==1
        counts_1(:,3) = hist(d_ind,1:nbins_max);
    elseif K==2
        counts_2(:,3) = hist(d_ind,1:nbins_max);
    elseif K==4
        counts_4(:,3) = hist(d_ind,1:nbins_max);
    end

end
%%

counts_1
counts_2
counts_4

