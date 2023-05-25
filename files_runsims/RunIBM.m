%% Test: run the first few simulations

clear all; clc; close all
% written by Allison K. Shaw (ashaw@umn.edu)
% updated: June 2018

% NSLOTS = str2num(getenv('NSLOTS'));
% maxNumCompThreads(NSLOTS);

clear all

saveto_filepath = '../output_simulations/20230524_test';

%BIOLOGICAL PARAMETERS
gflag = 0;     % whether (1) or not (0) to show plots during simulation
b = 10;        % offspring produced per individual
del = 0.001;   % fraction of dispersal probability to move during mutation
p = .75;      % probability of surviving dispersal

%ENVIRONMENT PARAMETERS
nbins_env = 12; % max number of dispersal bins used when creating env (>=2)
nbins = 6; % how many dispersal bins to actually use (2<=nbins_use<=nbins_max)
eflag = 2; % which environment to use: 1=unbounded, 2=bounded, 3=reef
S = 32^2;      % number of sites in the environment
sx = 2; % number of sites in the x-dimension of the environment
sy = S/sx;
    
%BIOLOGICAL PARAMETERS CHANGED
G = 5; % number of total generations to simulate
nmax = 0;  % maximum larval navigation distance (behavior)

%INITIAL CONDITIONS
    % option 1: specify starting displacement kernel
    v = [0.5 0.25 0.2 0.04 0.01 0]; % starting displacement kernel
    pop_init = fn_create_initial_pop(v,eflag,sx,sy,nbins_env,nbins,nmax,b);
    
    % option 2: default, uniform displacement kernel
    % pop_init=0;
    
    % option 3: specify starting population
    % (e.g., from end of a previous simulation)
    % load('../output_simulations/20220805/pop_nmax=0.mat');
    % pop_init=pop;
    % clear pop

fn_IBM_dispersal(gflag,eflag,sx,sy,nbins_env,nbins,nmax,G,b,p,del,pop_init,saveto_filepath)
