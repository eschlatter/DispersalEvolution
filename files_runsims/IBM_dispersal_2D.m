function [] = IBM_dispersal_2D(gflag,eflag,sx,sy,nbins_env,nbins,nmax,G,b,p,del)
% written by Allison K. Shaw (ashaw@umn.edu)
% updated: June 2018
%
% IBM version of Hamilton & May model with spatially-explicit patches in
% 2-D
% Two steps to the dispersal process: individuals go from natal patch to
%   landing patch (using dmap), then move from landing patch to final
%   patch (using tmap)
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

saveto_filepath = '../output_simulations/20220516';

if nbins < 2; error('nbins_use must be at least 2'); end
if nbins > nbins_env; error('nbins_env must be bigger than nbins'); end

%-----LOAD-ENVIRONMENT----------------------------------------------------%
    if eflag==1
        load(strcat(['../output_environments/env_unbounded_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins_env) '.mat']))
    elseif eflag==2
        load(strcat(['../test_output_environments/env_bounded_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins_env) '_nmax=' num2str(nmax) '.mat']))
    elseif eflag==3
        load(strcat(['../output_environments/env_reef_nbins=' num2str(nbins_env) '_nmax=' num2str(nmax) '_bmin=' num2str(b(1)) '_bmax=' num2str(b(2)) '.mat']))
    elseif eflag==4
        load(strcat(['../output_environments/env_bounded_het_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins_env) '_nmax=' num2str(nmax) '_bmin=' num2str(b(1)) '_bmax=' num2str(b(2)) '.mat']))
    end
%-----LOAD-ENVIRONMENT----------------------------------------------------%


%-----INITALIZATION-------------------------------------------------------%
    N0 = round(K*S); % initial number of individuals - at capacity

    % create matrix to hold information for all individuals
    pop = zeros(N0,nbins+1);
    pop(:,1:(nmax+2)) = (1/(nmax+2))*ones(N0,nmax+2); % assign starting dispersal bin values (set equal prob of each distance up to navigation distance + 1)
    %pop(:,1) = ones(N0,1); % assign starting dispersal bin values (set prob stay = 1, others zero)
    pop(:,nbins+1) = repmat(1:S,1,K); % assign evenly to starting patches

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
    N = size(pop,1);  % count number of individuals present
    dtime_avg(g,:) = mean(pop(:,1:nbins),1);   % average dispersal params
    dtime_std(g,:) = std(pop(:,1:nbins),[],1); % std. of dispersal params

    %-----REPRODUCTION-----%
        % each individual produces b offspring with same patch and dispersal
        % parameters as parent
        off = [];
        for j = 1:size(pop,1)
            off = [off; repmat(pop(j,:),[bvec(pop(j,nbins+1)),1])];
        end
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

    pat_disp = zeros(Noff,1); % to hold the patch where each offspring lands after dispersal (displacement)

    ind = find(d_ind==1); % which offspring did not disperse
    pat_disp(ind) = dmap{1}(off(ind,nbins+1)); % map natal patch from dmap
    clear ind

    % for each of the possible dispersal distances
    for j = 2:nbins

        ind = find(d_ind==j); % find all offspring that drew j-1 dipersal distance

        for i = ind'
            % for offspring that survived dispersal
            if surv(i)==1
                x = dmap{j}(off(i,nbins+1),:); % possible patches to disperse to
                y = randi(length(x)); % pick one at random
                pat_disp(i) = x(y); % save landing patch
            end
        end

    end

    % remove offspring who died during dispersal here
    off(pat_disp==0,:) = []; % remove offspring who died during dispersal
    pat_disp(pat_disp==0) = [];
    fitness(g,2) = length(pat_disp); % record number of offspring left after dispersal mortality

    pat_sett = tmap(pat_disp); % where offspring settle after dispersal

    off(:,nbins+1) = pat_sett; % save new locations

    % find offspring who didn't settle anywhere and remove
    off(off(:,nbins+1)==0,:) = []; % remove offspring who died during dispersal
    fitness(g,3) = length(off); % record number of offspring left after navigation



    clear srand surv drand d_ind ind i j pat_disp
    %-----MUTATION-AND-DISPERSAL-----%


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
    end
    clear fullind patch_ind i
    % kill off all adults and just save offspring as new population
    pop = off;
    %-----COMPETITION-----%

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

%-----PLOT-RESULTS--------------------------------------------------------%
    figure(2);clf
    % plot the evolved dispersal bin values across generations
    %subplot(2,2,3)
    hold on
    errorbar(dtime_avg,dtime_std)
    axis([0 G 0 1.1])
    xlabel('generation number')
    ylabel('probability')
    title('dispersal bin values')
    legend(string(1:nbins))
    % save figure
    if eflag==1
        saveas(2,strcat([saveto_filepath '/IBMfig_unbounded_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins) '_nmax=' num2str(nmax) '_del=' num2str(del) '_b=' num2str(b) '_p=' num2str(p) '.jpg']))
    elseif eflag==2
        saveas(2,strcat([saveto_filepath '/IBMfig_bounded_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins) '_nmax=' num2str(nmax) '_del=' num2str(del) '_b=' num2str(b) '_p=' num2str(p) '.jpg']))
    elseif eflag==3
        saveas(2,strcat([saveto_filepath '/IBMfig_reef_nbins=' num2str(nbins) '_nmax=' num2str(nmax) '_bmin=' num2str(bmin) '_bmax=' num2str(bmax) '_del=' num2str(del) '_p=' num2str(p) '.jpg']))
    elseif eflag==4
        saveas(2,strcat([saveto_filepath '/IBMfig_bounded_het_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins_env) '_nmax=' num2str(nmax) '_bmin=' num2str(bmin) '_bmax=' num2str(bmax) '_del=' num2str(del) '_p=' num2str(p) '.jpg']))
    end

    % plot the kernel of population mean displacement probabilities
    figure(3);clf
    kern = sum(pop(:,1:nbins),1)/length(pop);
    bar(0:(nbins-1),kern,'k')
    axis([-0.5 29.5 0 1])
    ylabel('Population mean probability')
    xlabel('Distance')
    title(sprintf('Dispersal Kernel, nmax = %g',nmax))

    writematrix(kern,sprintf('%s/kernel_WithNav_nbins=30_nmax = %g.csv',saveto_filepath,nmax))
    saveas(3, sprintf('%s/fig_WithNav_nbins=30_nmax=%g.jpg',saveto_filepath,nmax))
%-----PLOT-RESULTS--------------------------------------------------------%
