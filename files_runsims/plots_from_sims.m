% look at simulation results
clear

load('C:\Users\eschlatter\Dropbox\DispersalEvolution\output_simulations\20220513\IBM_bounded_sx=2_sy=512_nbins=30_nmax=2_del=0.001_b=10_p=1.mat')

% fitness over time
% (fitness = fraction of larvae that survive after dispersal and navigation)
plot(1:G,fitness(:,3)./fitness(:,1))
xlabel('Generation')
ylabel('Fraction of larve surviving dispersal and navigation')

