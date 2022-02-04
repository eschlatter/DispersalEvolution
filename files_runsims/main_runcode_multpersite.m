clear all; clc; close all
% written by Allison K. Shaw (ashaw@umn.edu)
% updated: June 2018

gflag = 0;     % whether (1) or not (0) to show plots during simulation
del = 0.001;   % fraction of dispersal probability to move during mutation
nbins_env = 31; % max number of dispersal bins used when creating env (>=2)
p = 1; % probability of surviving dispersal
sx = []; sy = []; % don't need to specify these for reef env
nbins = 31; % how many dispersal bins to actually use (2<=nbins_use<=nbins_max)
nmax = 0; % maximum larval navigation distance (behavior)
G = 200000; % number of total generations to simulate
b0 = 10;        % offspring produced per individual
eflag = 3; %which environment to use: 1=unbounded, 2=bounded, 3=reef, 4=heterogeneous

% seascape with no patchiness
bmin = b0; % minimum offspring per site in the environment
bmax = b0; % maximum offspring per site in the environment
b = [bmin bmax];

K = 2; % carrying capacity per patch
IBM_dispersal_2D_wK(gflag,eflag,sx,sy,nbins_env,nbins,nmax,G,b,p,del,K)

K = 4; % carrying capacity per patch
IBM_dispersal_2D_wK(gflag,eflag,sx,sy,nbins_env,nbins,nmax,G,b,p,del,K)
