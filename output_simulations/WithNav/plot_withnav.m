clear

S = 32^2;      % number of sites in the environment    
del = 0.001;   % fraction of dispersal probability to move during mutation
b = 10;        % offspring produced per individual
pvec1 = 1; % probability of surviving dispersal
nbins_max = 30; % how many dispersal bins to actually use (2<=nbins_use<=nbins_max)
sx = 2;
sy = S/sx;

for ii = 1:4

    nmax = ii-1;  % maximum larval navigation distance (behavior)

    load(strcat(['IBM_bounded_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins_max) '_nmax=' num2str(nmax) '_del=' num2str(del) '_b=' num2str(b) '_p=' num2str(pvec1) '.mat']))
    
    figure(ii)
%     tiledlayout(1,2)
%     
%     nexttile
%     
%     errorbar(dtime_avg,dtime_std)
%     axis([0 G 0 1.1])
%     xlabel('Generation')
%     ylabel('Probability')
%     title('Dispersal Evolution')
%     legend(string((1:nbins)-1))
%     
%     nexttile
%     
    kern = sum(pop(:,1:30),1)/length(pop);
    bar(0:29,kern,'k')
    axis([-0.5 29.5 0 1])
    ylabel('Population mean probability')
    xlabel('Distance')
    title(sprintf('Dispersal Kernel, nmax = %g',nmax))

    writematrix(kern,sprintf('kernel_nmax = %g.csv',nmax))
    %saveas(ii, sprintf('fig_WithNav_nbins=30_nmax=%g.jpg',nmax))
end