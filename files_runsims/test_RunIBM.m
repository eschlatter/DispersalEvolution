%% Test: run the first few simulations

clear all; clc; close all
% written by Allison K. Shaw (ashaw@umn.edu)
% updated: June 2018

% NSLOTS = str2num(getenv('NSLOTS'));
% maxNumCompThreads(NSLOTS);

clear all

%BIOLOGICAL PARAMETERS
gflag = 0;     % whether (1) or not (0) to show plots during simulation
b = 10;        % offspring produced per individual
del = 0.001;   % fraction of dispersal probability to move during mutation
p = 1;      % probability of surviving dispersal

%ENVIRONMENT PARAMETERS
nbins_env = 12; % max number of dispersal bins used when creating env (>=2)
nbins = 12; % how many dispersal bins to actually use (2<=nbins_use<=nbins_max)
eflag = 2; % which environment to use: 1=unbounded, 2=bounded, 3=reef
S = 32^2;      % number of sites in the environment
sx = 2; % number of sites in the x-dimension of the environment
sy = S/sx;
    
G = 5; % number of total generations to simulate

% larval navigation distance of 0, 1, 2 or 3
nmax = 0;  % maximum larval recruitment distance (behavior)

new_IBM_dispersal_2D(gflag,eflag,sx,sy,nbins_env,nbins,nmax,G,b,p,del)
