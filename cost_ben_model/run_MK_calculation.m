clear

%KERNEL PARAMETERS
nbins = 6; % how many dispersal bins to use (2<=nbins)
ndiv = 6; % increments by which to shift when generating kernels

%BIOLOGICAL PARAMETERS
b = 10;        % offspring produced per individual
p = 1;      % probability of surviving dispersal

%ENVIRONMENT PARAMETERS
eflag = 2; % which environment to use: 1=unbounded, 2=bounded, 3=reef
S = 32^2;      % number of sites in the environment
sx = 2; % number of sites in the x-dimension of the environment
sy = S/sx;
    
% larval navigation distance of 0, 1, 2 or 3
nmax = 0;  % maximum larval recruitment distance (behavior)

load('weibull_kerns.mat')
saveto_filepath = 'weibull_MK';

fn_calculate_MK(kerns,b,p,nbins,eflag,sx,sy,nmax,saveto_filepath)