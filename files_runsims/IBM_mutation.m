function [] = IBM_mutation(S,nbins,G,del)

% INPUTS:
%   S = number of sites
%   nbins = number of possible displacement distances
%   G = number of generations
%   del = amount of mutation

rng('shuffle') % seed the random number generator from computer clock

K = 1;         % carrying capacity per patch
nmax = 0;       % navigation ability

saveto_filepath = '../output_simulations/20221021_test';

if nbins < 2; error('nbins_use must be at least 2'); end

%-----INITALIZATION-------------------------------------------------------%
    N0 = round(K*S); % initial number of individuals - at capacity

    % create matrix to hold information for all individuals
    pop = zeros(N0,nbins+1);
    pop(:,1:(nmax+2)) = (1/(nmax+2))*ones(N0,nmax+2); % assign starting dispersal bin values (set equal prob of each distance up to navigation distance + 1)
    
    % matrices to record population parameters in, over time
    dtime_avg = zeros(G,nbins); % population mean of dispersal parameters
    dtime_std = zeros(G,nbins); % pop. standard deviation of dispersal parammeters

%-----INITALIZATION-------------------------------------------------------%

tic

%-----SIMULATE------------------------------------------------------------%
g=0; % this counts the number of generations that have passed

while g<G && size(pop,1)>0 % loop over generations (only while population not extinct)
    g = g+1 % step generation forward

    % record statistics
    dtime_avg(g,:) = mean(pop(:,1:nbins),1);   % average dispersal params
    dtime_std(g,:) = std(pop(:,1:nbins),[],1); % std. of dispersal params

    %-------------------------MUTATION--------------------%

    for j = 1:N0

        % mutate dispersal strategy slightly by taking or adding del from dispersal bin
        % if del amount isn't remaining in dispersal bin, take everything
        % (i.e. bound dispersal probability at zero)
        ind = randperm(nbins,2);         % select two bins at random
        val = min(pop(j,ind(1)),del);          % take del of prob, if available
        pop(j,ind) = pop(j,ind) + [-val +val]; % move del from first to second

    end

    clear ind j val
    %----------------------------MUTATION-----------------%
end

    % save output as mat file
    save(strcat([saveto_filepath '/IBM_mutation_S=' num2str(S) '_nbins=' num2str(nbins) '_G=' num2str(G) '_del=' num2str(del) '.mat']))

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
    saveas(2,strcat([saveto_filepath '/IBMfig_mutation_S=' num2str(S) '_nbins=' num2str(nbins) '_G=' num2str(G) '_del=' num2str(del) '.jpg']))

    % plot the kernel of population mean displacement probabilities
    figure(3);clf
    kern = sum(pop(:,1:nbins),1)/length(pop);
    bar(0:(nbins-1),kern,'k')
    axis([-0.5 nbins-0.5 0 1])
    ylabel('Population mean probability')
    xlabel('Distance')
    title(sprintf('Dispersal Kernel, nmax = %g',nmax))

    writematrix(kern,sprintf('%s/kernel_mutation_S=%g_nbins=%g_G=%g_del=%g.csv',saveto_filepath,S,nbins,G,del))
    saveas(3, sprintf('%s/fig_kernel_mutation_S=%g_nbins=%g_G=%g_del=%g.jpg',saveto_filepath,S,nbins,G,del))
%-----PLOT-RESULTS--------------------------------------------------------%
end
