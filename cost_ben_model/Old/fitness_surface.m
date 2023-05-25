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
nmax = 1;  % maximum larval recruitment distance (behavior)

% --------------- Generate kernels ------------------ %

bars = nchoosek(1:(nbins+ndiv-1),(nbins-1)); % locations of bars in "bars and stars" proof
bars = [zeros(length(bars),1) bars (nbins+ndiv)*ones(length(bars),1)];

kernels = ones(length(bars),nbins);

for j=1:length(bars)
    for i=1:nbins
        kernels(j,i) = bars(j,i+1)-bars(j,i)-1;
    end
end

kernels = (1/ndiv)*kernels;

% --------------- Get costs for each kernel ------------------- %
%outputs = array2table(zeros(length(kernels),3),'VariableNames',{'Fitness','Kin Comp','Mortality'});
out_fitness = zeros(length(kernels),1);
out_kincost = zeros(length(kernels),1);
out_mortcost = zeros(length(kernels),1);
dist_fitness = zeros(length(kernels),nbins);
dist_kincost = zeros(length(kernels),nbins);
dist_mortcost = zeros(length(kernels),nbins);

tic
for i=1:size(kernels,1)
    v = kernels(i,:);
    [M,K]=fn_one_costben_sim(b,p,nbins,eflag,sx,sy,nmax,v);
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

save('nmax_1/fitness_surface.mat',"out_fitness","out_kincost","out_mortcost","dist_fitness","dist_kincost","dist_mortcost","kernels")


