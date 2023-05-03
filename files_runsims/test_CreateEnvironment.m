clear all; clc; close all
% written by Allison K. Shaw (ashaw@umn.edu)
% updated by E Schlatter, April 2022

NSLOTS = str2num(getenv('NSLOTS'));
maxNumCompThreads(NSLOTS);

% STEP 1: create environments, before moving on to IBMs

    S = 32^2;      % number of sites in the environment
    b = 10;        % offspring produced per individual
    sx = 2; % number of sites in the x-dimension of the environment
    sy = S/sx;
    create_env_bounded(sx,sy,b)

%% BOUNDED ENVIRONMENTS

    nbins_env = 10; % max number of dispersal bins the env can support (must be >= 2)
    S = 32^2;      % number of sites in the environment
    b = 10;        % offspring produced per individual

    % bounded 2x512 environment WITHOUT larval navigation (nmax = 0)
    sx = 2; % number of sites in the x-dimension of the environment
    sy = S/sx;
    nmax = 2;      % maximum larval navigation distance (behavior)
    create_env_bounded(sx,sy,nbins_env,nmax,b)
    close all; clear sx sy nmax
    
    % bounded 2x512 environment WITH larval navigation (nmax = 1)
    sx = 2; % number of sites in the x-dimension of the environment
    sy = S/sx;
    nmax = 1;      % maximum larval navigation distance (behavior)
    create_env_bounded(sx,sy,nbins_env,nmax,b)
    close all; clear sx sy nmax

    % bounded 2x512 environment WITH larval navigation (nmax = 2)
    sx = 2; % number of sites in the x-dimension of the environment
    sy = S/sx;
    nmax = 2;      % maximum larval navigation distance (behavior)
    create_env_bounded(sx,sy,nbins_env,nmax,b)
    close all; clear sx sy nmax

        % bounded 2x512 environment WITH larval navigation (nmax = 3)
    sx = 2; % number of sites in the x-dimension of the environment
    sy = S/sx;
    nmax = 3;      % maximum larval navigation distance (behavior)
    create_env_bounded(sx,sy,nbins_env,nmax,b)
    close all; clear sx sy nmax
    
    clear S
