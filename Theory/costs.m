load('../output_simulations/20220805/IBM_bounded_sx=2_sy=512_nbins=12_nmax=0_del=0.001_b=10_p=1.mat')
%% basic costs figure

% vectors of parameters

%v = 0.1*ones(1,10); % displacement strategy
v = mean(pop(:,1:12));
%v = [0.25 0.25 0.25 0.25 0 0 0 0 0 0 0 0];
D = length(v); % number of displacement bins
n = [1 4*(1:D-1)]; % number of sites at each distance
g = [1 3 4*ones(1,D-2)]; % number of viable sites at each distance
c = zeros(1,D); % displacement process mortality

% calculate direct and indirect costs of displacement to each distance

% direct costs:
% (risk to larva of mortality due to i) displacement process or ii) landing
% in a nonviable site)
m = c+(1-c).*(1-g./n);

% indirect costs
% (due to kin competition among full-sibs)
S = sum(g.*v.*(1-c)./n); %denominator of the quantity
k = (1/S)*v.*(1-c)./n;

% total benefit
% product of direct and indirect benefits (one minus each cost)
prod = (1-m).*(1-k); 

figure
hold on
bar(v,'w')
plot(1:D,m,'.-',1:D,k,'.-',1:D,prod,'.-')
legend('kernel','direct cost','indirect cost','total benefit')
xlabel('Distance')
hold off

%% error bars

%get max and min for last ### generations
end_seq = dtime_avg(70001:100000,:);
[~,max_ind] = max(end_seq(:,1));
max_kern = end_seq(max_ind,:); %kernel with the max value in bin 0
[~,min_ind] = min(end_seq(:,1));
min_kern = end_seq(min_ind,:); %kernel with the min value in bin 0

%for max

max_out = kernelcosts(max_kern,1);
min_out = kernelcosts(min_kern,1);
mean_ind_cost = mean([max_out(2,:);min_out(2,:)]);
mean_benefit = mean([max_out(3,:);min_out(3,:)]);
ind_cost_error = abs(max_out(2,:)-min_out(2,:));
benefit_error = abs(max_out(3,:)-min_out(3,:));

D = length(max_kern);
n = [1 4*(1:D-1)]; % number of sites at each distance
g = [1 3 4*ones(1,D-2)]; % number of viable sites at each distance
c = zeros(1,D); % displacement process mortality
% direct costs:
% (risk to larva of mortality due to i) displacement process or ii) landing
% in a nonviable site)
m = c+(1-c).*(1-g./n);

figure
hold on
bar(max_kern,'FaceColor','none')
bar(min_kern,'FaceColor','none')
plot(1:D,m,'.-r')
errorbar(1:D,mean_ind_cost,ind_cost_error/2,'Color','#EDB120')
errorbar(1:D,mean_benefit,benefit_error/2,'Color','#7E2F8E')
legend('kernel1','kernel2','direct cost','indirect cost','total benefit')
xlabel('Distance')
hold off
