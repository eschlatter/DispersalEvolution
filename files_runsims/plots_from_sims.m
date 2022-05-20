% look at simulation results
clear

load('C:\Users\eschlatter\Dropbox\DispersalEvolution\output_simulations\20220516\IBM_bounded_sx=2_sy=512_nbins=12_nmax=2_del=0.001_b=10_p=1.mat')

totals = [sum(kernel_displacement,2) sum(kernel_dispersal,2) sum(kernel_recruitment,2)];




% fitness over time
% (fitness = fraction of larvae that survive after dispersal and navigation)
plot(1:G,fitness(:,3)./fitness(:,1))
xlabel('Generation')
ylabel('Fraction of larve surviving dispersal and navigation')

% plot the evolved dispersal bin values across generations
figure
hold on
errorbar(dtime_avg,dtime_std)
axis([0 G 0 1.1])
xlabel('generation number')
ylabel('probability')
title('dispersal bin values')
legend(string((1:nbins)-1))
hold off

% plot the evolving displacement kernel across generations
figure
plot(kernel_displacement)
legend(string((1:nbins)-1))
title('Displacement Kernel')

% plot the evolving dispersal kernel across generations
figure
plot(kernel_dispersal)
legend(string((1:nbins)-1))
title('Dispersal Kernel')

% plot the evolving recruitment kernel across generations
figure
plot(kernel_recruitment)
legend(string((1:nbins)-1))
title('Recruitment Kernel')