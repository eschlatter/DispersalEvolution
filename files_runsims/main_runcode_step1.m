clear all; clc; close all
% written by Allison K. Shaw (ashaw@umn.edu)
% updated: June 2018

% STEP 1: create all possible environments, before moving on to IBMs


%% WITHOUT BOUND ENVIRONMENT (1)

    nbins_env = 4; % max number of dispersal bins the env can support (must be >= 2)
    S = 32^2;      % number of sites in the environment
    b = 10;        % offspring produced per individual

    % unbounded 32x32 environment with no larval navigation
    sx = 32; % number of sites in the x-dimension of the environment
    sy = S/sx;
    create_env_unbounded(sx,sy,nbins_env,b)
    close all; clear sx sy nmax
    
%% WITH BOUND ENVIRONMENTS (3)

    nbins_env = 4; % max number of dispersal bins the env can support (must be >= 2)
    S = 32^2;      % number of sites in the environment
    b = 10;        % offspring produced per individual

    % bounded 32x32 environment with no larval navigation
    sx = 32; % number of sites in the x-dimension of the environment
    sy = S/sx;
    nmax = 0;      % maximum larval navigation distance (behavior)
    create_env_bounded(sx,sy,nbins_env,nmax,b)
    close all; clear sx sy nmax
    
    % bounded 8x128 environment with no larval navigation
    sx = 8; % number of sites in the x-dimension of the environment
    sy = S/sx;
    nmax = 0;      % maximum larval navigation distance (behavior)
    create_env_bounded(sx,sy,nbins_env,nmax,b)
    close all; clear sx sy nmax
    
    % bounded 2x512 environment with no larval navigation
    sx = 2; % number of sites in the x-dimension of the environment
    sy = S/sx;
    nmax = 0;      % maximum larval navigation distance (behavior)
    create_env_bounded(sx,sy,nbins_env,nmax,b)
    close all; clear sx sy nmax
    
    clear S


%% WITH BOUND ENVIRONMENTS, HETEROGENEOUS (3)

% bounded 2x512 environment with no larval navigation
    nbins_env = 4; % max number of dispersal bins the env can support (must be >= 2)
    S = 32^2;      % number of sites in the environment
    sx = 2; % number of sites in the x-dimension of the environment
    sy = S/sx;
    nmax = 0;      % maximum larval navigation distance (behavior)
    b = 10;        % offspring produced per individual

    % low-level patchiness (between 2b/3 and 4b/3)
    bmin = round(2*b/3); % minimum offspring per site in the environment
    bmax = round(4*b/3); % maximum offspring per site in the environment
    create_env_heterogeneous(sx,sy,nbins_env,nmax,bmin,bmax)
    clear bmin bmax

    % mid-level patchiness (between b/3 and 5b/3)
    bmin = round(b/3); % minimum offspring per site in the environment
    bmax = round(5*b/3); % maximum offspring per site in the environment
    create_env_heterogeneous(sx,sy,nbins_env,nmax,bmin,bmax)
    clear bmin bmax

    % high-level patchiness (between 0 and 2b)
    bmin = round(0); % minimum offspring per site in the environment
    bmax = round(2*b); % maximum offspring per site in the environment
    create_env_heterogeneous(sx,sy,nbins_env,nmax,bmin,bmax)
    clear bmin bmax


%% SEASCAPE, HETEROGENEOUS (3)
    nbins_env = 31; % max number of dispersal bins the env can support (must be >= 2)
    nmax = 0;      % maximum larval navigation distance (behavior)
    b = 10;        % offspring produced per individual
    
    % seascape with no patchiness
    bmin = b; % minimum offspring per site in the environment
    bmax = b; % maximum offspring per site in the environment
    create_env_reef(nbins_env,nmax,bmin,bmax)
    clear bmin bmax

    
    % seascape with low-level patchiness (between 2b/3 and 4b/3)
    bmin = round(2*b/3); % minimum offspring per site in the environment
    bmax = round(4*b/3); % maximum offspring per site in the environment
    create_env_reef(nbins_env,nmax,bmin,bmax)
    clear bmin bmax

    % seascape with mid-level patchiness (between b/3 and 5b/3)
    bmin = round(b/3); % minimum offspring per site in the environment
    bmax = round(5*b/3); % maximum offspring per site in the environment
    create_env_reef(nbins_env,nmax,bmin,bmax)
    clear bmin bmax

    % seascape with high-level patchiness (between 0 and 2b)
    bmin = round(0); % minimum offspring per site in the environment
    bmax = round(2*b); % maximum offspring per site in the environment
    create_env_reef(nbins_env,nmax,bmin,bmax)
    clear bmin bmax

