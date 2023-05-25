clear

saveto_filepath = 'weibull_MK';

%KERNEL PARAMETERS
nbins = 6; % how many dispersal bins to use (2<=nbins)
ndiv = 6; % increments by which to shift when generating kernels

%BIOLOGICAL PARAMETERS
b = 10;        % offspring produced per individual
p = 1;      % probability of surviving dispersal
nmax = 0;  % maximum larval recruitment distance (behavior)

%ENVIRONMENT PARAMETERS
eflag = 2; % which environment to use: 1=unbounded, 2=bounded, 3=reef
S = 32^2;      % number of sites in the environment
sx = 2; % number of sites in the x-dimension of the environment
sy = S/sx;
    

%KERNEL
    %-----Parametric version (load a set of kernels)
    load('weibull_kerns.mat')
    
%     %-----Nonparametric version (calculate from nbins and ndiv)
%     bars = nchoosek(1:(nbins+ndiv-1),(nbins-1)); % locations of bars in "bars and stars" proof
%     bars = [zeros(length(bars),1) bars (nbins+ndiv)*ones(length(bars),1)];
%     
%     kerns = ones(length(bars),nbins);
%     
%     for j=1:length(bars)
%         for i=1:nbins
%             kerns(j,i) = bars(j,i+1)-bars(j,i)-1;
%         end
%     end
%     
%     kerns = (1/ndiv)*kerns;
%     %-----Nonparametric version

fn_many_costben_sims(kerns,b,p,nbins,eflag,sx,sy,nmax,saveto_filepath)