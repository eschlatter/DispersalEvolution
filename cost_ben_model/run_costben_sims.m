clear

%BIOLOGICAL PARAMETERS
b = 10;        % offspring produced per individual
p = 1;      % probability of surviving dispersal

%ENVIRONMENT PARAMETERS
nbins = 6; % how many dispersal bins to use (2<=nbins)
eflag = 2; % which environment to use: 1=unbounded, 2=bounded, 3=reef
S = 32^2;      % number of sites in the environment
sx = 2; % number of sites in the x-dimension of the environment
sy = S/sx;
    
% larval navigation distance of 0, 1, 2 or 3
nmax = 0;  % maximum larval recruitment distance (behavior)

v = [0.5 0.25 0.2 0.04 0.01 0];

%-----------Run the function----------------%
[M,K,K2]=fn_costben_sim(b,p,nbins,eflag,sx,sy,nmax,v);
B = (1-M).*(1-K2); % total benefit


%-------------Plot output---------------------%
figure
hold on
bar(v,'w')
plot(1:nbins,M,'.-','Color','#D95319','LineWidth',1,'MarkerSize',10)
plot(1:nbins,K,'.-','Color','#EDB120','LineWidth',1,'MarkerSize',10)
plot(1:nbins,K2,'--','Color','#EDB120','LineWidth',1,'MarkerSize',10)
plot(1:nbins,B,'--','Color','#7E2F8E','LineWidth',1,'MarkerSize',10)
legend('kernel','mortality','kin competition','total benefit')
xlabel('Distance')
xticks(1:nbins)
xticklabels(string(0:nbins-1))
ylim([0 1])
%title(sprintf('nmax=%g',nmax))
legend('off')
hold off

%------------------sum over all distances-------------------%
M(isnan(M))=0; % NaN entries mean no larvae displaced that distance,
K(isnan(K))=0; % so no contribution to costs/benefits
K2(isnan(K2))=0;
BT = v*B;
KT = v*K;
K2T = v*K2;
MT = v*M;