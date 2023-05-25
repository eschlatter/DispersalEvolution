function []=fn_calculate_MK(kernels,b,p,nbins,eflag,sx,sy,nmax,saveto_filepath)
% takes a matrix of kernels (each row a kernel, each column a distance
% probability) plus simulation parameters, and outputs mortality and kin
% competition costs and fitness for each kernel

out_fitness = zeros(size(kernels,1),1);
out_kincost = zeros(size(kernels,1),1);
out_mortcost = zeros(size(kernels,1),1);
dist_fitness = zeros(size(kernels,1),nbins);
dist_kincost = zeros(size(kernels,1),nbins);
dist_mortcost = zeros(size(kernels,1),nbins);

tic
for i=1:size(kernels,1)
    i
    v = kernels(i,:);
    [M,K]=fn_costben_sim(b,p,nbins,eflag,sx,sy,nmax,v);
    B = (1-M).*(1-K); % total benefit
    dist_fitness(i,:) = B;
    B(isnan(B))=0; % NaN entries mean no larvae displaced that distance,
    out_fitness(i) = v*B;
    dist_kincost(i,:) = K;
    dist_mortcost(i,:) = M;
    M(isnan(M))=0; % NaN entries mean no larvae displaced that distance,
    K(isnan(K))=0; % so no contribution to costs/benefits
    out_kincost(i) = v*K;
    out_mortcost(i) = v*M;
    clear K M B v
end
toc

save(strcat([saveto_filepath,'.mat']),"out_fitness","out_kincost","out_mortcost","dist_fitness","dist_kincost","dist_mortcost","kernels")