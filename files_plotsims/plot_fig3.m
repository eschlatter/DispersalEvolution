clear all; clc; close all
% written by Allison K. Shaw (ashaw@umn.edu)
% updated: October 2015
    
% PARAMETERS
    S = 32^2;      % number of sites in the environment
    nmax = 0;  % maximum larval navigation distance (behavior)
    del = 0.001;   % fraction of dispersal probability to move during mutation
    b = 10;        % offspring produced per individual
    pvec1 = 0.1:0.1:1; % probability of surviving dispersal
    nbins_max = 4; % how many dispersal bins to actually use (2<=nbins_use<=nbins_max)

figure(1); clf

subplot(3,1,1)
    sx = 32;
    sy = S/sx;
    sim_bounded = NaN(length(pvec1),nbins_max);
    for ii = 1:length(pvec1)
        load(strcat(['../output_simulations/IBM_bounded_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins_max) '_nmax=' num2str(nmax) '_del=' num2str(del) '_b=' num2str(b) '_p=' num2str(pvec1(ii)) '.mat']))
        sim_bounded(ii,:) = dtime_avg(end,:);
    end
    hold on
    bar(pvec1,sim_bounded,'stacked')
    axis([0 1.1 0 1])
    box on
    clear sim_bounded

subplot(3,1,2)
    sx = 8; % number of sites in the x-dimension of the environment
    sy = S/sx;
    sim_bounded = NaN(length(pvec1),nbins_max);
    for ii = 1:length(pvec1)
        load(strcat(['../output_simulations/IBM_bounded_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins_max) '_nmax=' num2str(nmax) '_del=' num2str(del) '_b=' num2str(b) '_p=' num2str(pvec1(ii)) '.mat']))
        sim_bounded(ii,:) = dtime_avg(end,:);
    end
    hold on
    bar(pvec1,sim_bounded,'stacked')
    axis([0 1.1 0 1])
    box on    
    clear sim_bounded

subplot(3,1,3)
    sx = 2; % number of sites in the x-dimension of the environment
    sy = S/sx;
    sim_bounded = NaN(length(pvec1),nbins_max);
    for ii = 1:length(pvec1)
        load(strcat(['../output_simulations/IBM_bounded_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins_max) '_nmax=' num2str(nmax) '_del=' num2str(del) '_b=' num2str(b) '_p=' num2str(pvec1(ii)) '.mat']))
        sim_bounded(ii,:) = dtime_avg(end,:);
    end
    hold on
    bar(pvec1,sim_bounded,'stacked')
    axis([0 1.1 0 1])
    box on
    xlabel('Survival probability')
