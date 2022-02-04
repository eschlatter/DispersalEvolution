clear all; clc; close all
% written by Allison K. Shaw (ashaw@umn.edu)
% updated: June 2018

% STEP 2: run IBMs for many different p values

gflag = 0;     % whether (1) or not (0) to show plots during simulation
b = 10;        % offspring produced per individual
del = 0.001;   % fraction of dispersal probability to move during mutation
pvec = 0.1:0.1:1; % probability of surviving dispersal

%% WITHOUT BOUND SIMULATIONS (3 x pvec)
    nbins_env = 4; % max number of dispersal bins used when creating env (>=2)

    % unbounded 32x32 environment with no larval navigation
    eflag = 1; %which environment to use: 1=unbounded, 2=bounded, 3=reef
    S = 32^2;      % number of sites in the environment
    sx = 32; % number of sites in the x-dimension of the environment
    sy = S/sx;
    nmax = 0;  % maximum larval recruitment distance (behavior)
    G = 100000; % number of total generations to simulate

    % 2 dispersal bins
    nbins = 2; % how many dispersal bins to actually use (2<=nbins_use<=nbins_max)
    for ii = 1:length(pvec)
        IBM_dispersal_2D(gflag,eflag,sx,sy,nbins_env,nbins,nmax,G,b,pvec(ii),del)
    end
    % 3 dispersal bins
    nbins = 3; % how many dispersal bins to actually use (2<=nbins_use<=nbins_max)
    for ii = 1:length(pvec)
        IBM_dispersal_2D(gflag,eflag,sx,sy,nbins_env,nbins,nmax,G,b,pvec(ii),del)
    end

    % 4 dispersal bins
    nbins = 4; % how many dispersal bins to actually use (2<=nbins_use<=nbins_max)
    for ii = 1:length(pvec)
        IBM_dispersal_2D(gflag,eflag,sx,sy,nbins_env,nbins,nmax,G,b,pvec(ii),del)
    end


%% BOUNDED SIMULATIONS (3 x pvec)
    nbins_env = 4; % max number of dispersal bins used when creating env (>=2)

    % bounded environment with no larval navigation, 4 dispersal bins
    eflag = 2; %which environment to use: 1=unbounded, 2=bounded, 3=reef
    S = 32^2;      % number of sites in the environment
    nmax = 0;  % maximum larval recruitment distance (behavior)
    nbins = 4; % how many dispersal bins to actually use (2<=nbins_use<=nbins_max)
    G = 100000; % number of total generations to simulate

    % bounded 32x32 environment with no larval navigation
    sx = 32; % number of sites in the x-dimension of the environment
    sy = S/sx;
    for ii = 1:length(pvec)
        IBM_dispersal_2D(gflag,eflag,sx,sy,nbins_env,nbins,nmax,G,b,pvec(ii),del)
    end
    
    % bounded 8x128 environment with no larval navigation
    sx = 8; % number of sites in the x-dimension of the environment
    sy = S/sx;
    for ii = 1:length(pvec)
        IBM_dispersal_2D(gflag,eflag,sx,sy,nbins_env,nbins,nmax,G,b,pvec(ii),del)
    end

    % bounded 2x512 environment with no larval navigation
    sx = 2; % number of sites in the x-dimension of the environment
    sy = S/sx;
    for ii = 1:length(pvec)
        IBM_dispersal_2D(gflag,eflag,sx,sy,nbins_env,nbins,nmax,G,b,pvec(ii),del)
    end



%%BOUNDED, HETEROGENEOUS (3 x pvec)
    nbins_env = 4; % max number of dispersal bins used when creating env (>=2)

    % bounded environment with no larval navigation, 4 dispersal bins
    eflag = 4; %which environment to use: 1=unbounded, 2=bounded, 3=reef, 4=heterogeneous
    S = 32^2;      % number of sites in the environment
    sx = 2; % number of sites in the x-dimension of the environment
    sy = S/sx;
    nmax = 0;  % maximum larval recruitment distance (behavior)
    nbins = 4; % how many dispersal bins to actually use (2<=nbins_use<=nbins_max)
    G = 100000; % number of total generations to simulate
    b0 = 10;        % offspring produced per individual

    % low-level patchiness (between 2b/3 and 4b/3)
    bmin = round(2*b0/3); % minimum offspring per site in the environment
    bmax = round(4*b0/3); % maximum offspring per site in the environment
    b = [bmin bmax];
    for ii = 1:length(pvec)
        IBM_dispersal_2D(gflag,eflag,sx,sy,nbins_env,nbins,nmax,G,b,pvec(ii),del)
    end

    % mid-level patchiness (between b/3 and 5b/3)
    bmin = round(b0/3); % minimum offspring per site in the environment
    bmax = round(5*b0/3); % maximum offspring per site in the environment
    b = [bmin bmax];
    for ii = 1:length(pvec)
        IBM_dispersal_2D(gflag,eflag,sx,sy,nbins_env,nbins,nmax,G,b,pvec(ii),del)
    end

    % high-level patchiness (between 0 and 2b)
    bmin = round(0); % minimum offspring per site in the environment
    bmax = round(2*b0); % maximum offspring per site in the environment
    b = [bmin bmax];
    for ii = 1:length(pvec)
        IBM_dispersal_2D(gflag,eflag,sx,sy,nbins_env,nbins,nmax,G,b,pvec(ii),del)
    end


%%SEASCAPE (4)
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
    IBM_dispersal_2D(gflag,eflag,sx,sy,nbins_env,nbins,nmax,G,b,p,del)

    % seascape with low-level patchiness (between 2b/3 and 4b/3)
    bmin = round(2*b0/3); % minimum offspring per site in the environment
    bmax = round(4*b0/3); % maximum offspring per site in the environment
    b = [bmin bmax];
    IBM_dispersal_2D(gflag,eflag,sx,sy,nbins_env,nbins,nmax,G,b,p,del)

    % seascape with mid-level patchiness (between b/3 and 5b/3)
    bmin = round(b0/3); % minimum offspring per site in the environment
    bmax = round(5*b0/3); % maximum offspring per site in the environment
    b = [bmin bmax];
    IBM_dispersal_2D(gflag,eflag,sx,sy,nbins_env,nbins,nmax,G,b,p,del)

    % seascape with high-level patchiness (between 0 and 2b)
    bmin = round(0); % minimum offspring per site in the environment
    bmax = round(2*b0); % maximum offspring per site in the environment
    b = [bmin bmax];
    IBM_dispersal_2D(gflag,eflag,sx,sy,nbins_env,nbins,nmax,G,b,p,del)


