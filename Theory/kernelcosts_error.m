function [dit, inds]=kernelcosts_error(kernel,dtime_avg,bin_to_min,nmax,showplot)
%calculates direct and indirect costs, and total benefit, of a given
%displacement strategy (for nmax=0,1,2)

v = kernel;
D = length(v); % number of displacement bins

c = zeros(1,D); % displacement process mortality
f = ones(1,D); % navigation survival

n = [1 4*(1:D-1)]; % number of sites at each distance

% number of viable sites at each distance (after navigation)
if nmax==0
    g = [1 3 4*ones(1,D-2)]; 
elseif nmax==1
    g = [1 4 7 8*ones(1,D-3)];
elseif nmax==2
    g = [1 4 8 11 12*ones(1,D-4)];
else
    sprintf('Not defined for this value of nmax')
end

% expected number of larvae displacing from one site to another, for sites distance i apart
L = 1./n.*(1-c).*v;

% Direct cost (mortality)
M = c+(1-c).*(1-g./n);

% Siblings at site at distance i (numerator of indirect cost)
if nmax==0
    S=L;
elseif nmax==1
    S=L+f(1).*[L(2:D) 0]; 
    %add the zero on the end because at the largest distance, no larvae can displace farther and then navigate back
elseif nmax==2
    Ssame = L+.25*f(1)*[L(2:D) 0]+.25*f(2)*(L+2*[L(3:D) 0 0])+f(2)*[L(3:D) 0 0];
    Sdiff = L+.25*f(1)*[L(2:D) 0]+.25*f(2)*(2*L+[L(3:D) 0 0])+f(2)*[L(3:D) 0 0];
    S = .5*(Ssame+Sdiff);
else
    sprintf('S not defined for this value of nmax')
end

% Total competitors at site at distance i from origin site (denominator of
% indirect cost)
if nmax==0
    T = sum(g.*L);
elseif nmax==1
    T = sum(g.*L)+f(1)*sum(g(1:D-1).*L(2:D));
elseif nmax==2
    T = sum(g.*L)+(.25*f(1)+.75*f(2))*sum(g(1:D-1).*L(2:D))+f(2)*sum(g(1:D-2).*L(3:D));
else
    sprintf('T not defined for this value of nmax')
end

% total indirect cost
K = S./T;

B = (1-M).*(1-K); % total benefit

dit = [M; K; B];

%get max and min for last ### generations
end_seq = dtime_avg(70001:100000,:);
[~,max_ind] = max(end_seq(:,bin_to_min));
max_kern = end_seq(max_ind,:); %kernel with the max value in bin 0
[~,min_ind] = min(end_seq(:,bin_to_min));
min_kern = end_seq(min_ind,:); %kernel with the min value in bin 0

inds=[max_ind min_ind]; %return indices of max and min used

max_out = kernelcosts(max_kern,nmax,0);
min_out = kernelcosts(min_kern,nmax,0);
mean_ind_cost = mean([max_out(2,:);min_out(2,:)]);
mean_benefit = mean([max_out(3,:);min_out(3,:)]);
error_ind_cost = abs(max_out(2,:)-min_out(2,:));
error_benefit = abs(max_out(3,:)-min_out(3,:));
  
if(showplot==1)

    hold on
    bar(max_kern,'FaceColor','none')
    bar(min_kern,'FaceColor','none')
    plot(1:D,M,'.-r')
    errorbar(1:D,mean_ind_cost,error_ind_cost/2,'Color','#EDB120')
    errorbar(1:D,mean_benefit,error_benefit/2,'Color','#7E2F8E')
    legend('kernel1','kernel2','direct cost','indirect cost','total benefit')
    xlabel('Distance')
    xticks(1:D)
    xticklabels(string(0:D-1))
    title(sprintf('nmax=%g',nmax))
    hold off

end

end