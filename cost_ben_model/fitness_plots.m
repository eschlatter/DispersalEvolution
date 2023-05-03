clear
load('fitness_surface.mat')
all = array2table([kernels out_mortcost out_kincost out_fitness],'VariableNames',{'v0','v1','v2','v3','v4','v5','Mortality','Kin Comp','Fitness'});

%% 2D: mort vs kin comp 

% basic scatterplot
plot(out_mortcost,out_kincost,'.')

% first, pull out just the points at the edge
edgepoints = [];
nslice = 20;
for i=1:nslice
    start = (i-1)/nslice;
    finish = i/nslice;
    slice = all((all.Mortality>=start & all.Mortality<finish),:);
    minslice = slice((slice.("Kin Comp")==min(slice.("Kin Comp"))),:);
    edgepoints = [edgepoints ; minslice];
    clear start finish slice minslice
end

clf
hold on
plot(out_mortcost,out_kincost,'.')
f = fit(edgepoints.Mortality,edgepoints.("Kin Comp"),"poly5");
plot(edgepoints.Mortality,edgepoints.("Kin Comp"),'o')
plot(f)
xlabel('Mortality Cost')
ylabel('Kin Competition Cost')
grid on
hold off



% add a trajectory
load('C:\Users\eschlatter\Dropbox\DispersalEvolution\output_simulations\20221103_1\dtime_1_nmax=0.mat')


%% 3D surface: mort and kin comp vs fitness

plot3(out_mortcost,out_kincost,out_fitness,'.') % just points
xlabel('Mortality Cost')
ylabel('Kin Competition Cost')
zlabel('Fitness')

f = fit([out_mortcost out_kincost],out_fitness,"poly11");
disp(f)
% fitness = 1 - mortcost - kincost
% aka, total cost = total mort cost (weighted avg, over all distances)
%                   + total kin cost (weighted avg over all distances)
plot(f,[out_mortcost out_kincost],out_fitness)
xlabel('Mortality Cost')
ylabel('Kin Competition Cost')
zlabel('Fitness')

% based on theory: fitness = 1 - mortcost - kincost + sum_d(v*m*k)
sum_d_vmk = sum(kernels .* dist_mortcost .* dist_kincost,2,'omitnan');
fitness_pred = sum_d_vmk + 1 - out_kincost - out_mortcost;
histogram(out_fitness-fitness_pred) %precise and unbiased -- the variation is due to stochasticity
fitness_pred2 = 1 - out_kincost - out_mortcost; %without that extra little term
histogram(out_fitness-fitness_pred2) %without that extra little term, fitness is systematically underestimated (by a small amount)

% what's the deal with that little term? How does it vary with kernel?
histogram(sum_d_vmk)
% what kernels have the largest? The ones with the mass concentrated into
% bins 1 and 2
biggest_terms = all(sum_d_vmk>0.045,:);
% what kernels have the smallest? The ones with most of their mass in bin
% 0, and the rest distributed among other bins
smallest_terms = all(sum_d_vmk<0.001,:);

rowss = kernels(:,1)==1/12;
plot3(out_mortcost(rowss,:),out_kincost(rowss,:),sum_d_vmk(rowss,:),'.') % just points

f = fit([kernels(:,1) out_kincost],sum_d_vmk,"cubicinterp");
plot(f,[kernels(:,1) out_kincost],sum_d_vmk)
xlabel('Mortality Cost')
ylabel('Kin Competition Cost')
zlabel('Little term')

% identify kernels in a region of the surface
stripe1 = all(out_kincost>0.8,:);
stripe2 = all(out_kincost<0.7 & out_kincost>0.5,:);
stripe3 = all(out_kincost<0.5 & out_kincost>0.3,:);

bestkernels = all(out_fitness>0.62,:);

% plot fitness vs kernel attributes
kernel_mean = kernels*(1:6)';

f = fit([kernel_mean out_kincost],out_fitness,"poly11");
plot(f,[kernel_mean out_kincost],out_fitness)
xlabel('Kernel Mean')
ylabel('Kin Competition Cost')
zlabel('Fitness')

