clear all; clc; close all
% written by Allison K. Shaw (ashaw@umn.edu)
% updated: October 2015
    
% PARAMETERS
    S = 32^2;      % number of sites in the environment
    sx = 32; % number of sites in the x-dimension of the environment
    sy = S/sx;
    nmax = 0;  % maximum larval navigation distance (behavior)
    del = 0.001;   % fraction of dispersal probability to move during mutation
    b = 10;        % offspring produced per individual
    pvec1 = 0.1:0.1:1; % probability of surviving dispersal
    pvec2 = 0:0.01:1; % probability of surviving dispersal

% ANALYIC
    % hamilton and may
    ana_hm = 1./(2-pvec2);
    % dispersing 0 or 1
    ana_unbound1 = 4./(8-3*pvec2);
    % dispersing 0, 1, or 2
    ana_unbound2(1,:) = 4./(24-11*pvec2); % distance 1
    ana_unbound2(2,:) = 2*ana_unbound2(1,:);       % distance 2
    % dispersing 0, 1, 2, or 3
    D = 3;
    ana_unbound3(1,:) = 4./(pvec2 + 2*D*(D+1)*(2-pvec2)); % distance 1
    ana_unbound3(2,:) = 2*ana_unbound3(1,:);              % distance 2
    ana_unbound3(3,:) = 3*ana_unbound3(1,:);              % distance 3

figure(1); clf

subplot(3,1,1)
    nbins_max = 2; %
    sim_unbounded = NaN(length(pvec1),nbins_max);
    for ii = 1:length(pvec1)
        load(strcat(['../output_simulations/IBM_unbounded_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins_max) '_nmax=' num2str(nmax) '_del=' num2str(del) '_b=' num2str(b) '_p=' num2str(pvec1(ii)) '.mat']))
       sim_unbounded(ii,:) = dtime_avg(end,:);
    end
    hold on
    bar(pvec1,[sim_unbounded zeros(10,2)],'stacked')
    plot(pvec2,1-ana_unbound1,'k-')
    plot(pvec2,1-ana_hm,'k:')
    axis([0 1.1 0 1])
    box on
    clear sim_unbounded

subplot(3,1,2)
    nbins_max = 3;
    sim_unbounded = NaN(length(pvec1),nbins_max);
    for ii = 1:length(pvec1)
        load(strcat(['../output_simulations/IBM_unbounded_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins_max) '_nmax=' num2str(nmax) '_del=' num2str(del) '_b=' num2str(b) '_p=' num2str(pvec1(ii)) '.mat']))
       sim_unbounded(ii,:) = dtime_avg(end,:);
    end
    hold on
    bar(pvec1,[sim_unbounded zeros(10,1)],'stacked')
    plot(pvec2,1-ana_unbound2(1,:)-ana_unbound2(2,:),'k-')
    plot(pvec2,1-ana_unbound2(2,:),'k--')
    plot(pvec2,1-ana_hm,'k:')
    axis([0 1.1 0 1])
    box on
    clear sim_unbounded

subplot(3,1,3)
    nbins_max = 4;
    sim_unbounded = NaN(length(pvec1),nbins_max);
    for ii = 1:length(pvec1)
        load(strcat(['../output_simulations/IBM_unbounded_sx=' num2str(sx) '_sy=' num2str(sy) '_nbins=' num2str(nbins_max) '_nmax=' num2str(nmax) '_del=' num2str(del) '_b=' num2str(b) '_p=' num2str(pvec1(ii)) '.mat']))
       sim_unbounded(ii,:) = dtime_avg(end,:);
    end
    hold on
    bar(pvec1,sim_unbounded,'stacked')
    plot(pvec2,1-ana_unbound3(1,:)-ana_unbound3(2,:)-ana_unbound3(3,:),'k-')
    plot(pvec2,1-ana_unbound3(2,:)-ana_unbound3(3,:),'k--')
    plot(pvec2,1-ana_unbound3(3,:),'k-.')
    plot(pvec2,1-ana_hm,'k:')
    box on
    xlabel('Survival probability')
    axis([0 1.1 0 1])
    
