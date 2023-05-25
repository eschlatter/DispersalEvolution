function comboplot2(dtime_avg,dtime_std,nmax,rep, intmin, intmax, endmin, endmax)
%plots kernel dynamics, plus costs figures for beginning, middle and end of
%evolutionary period

figure
clf
t=tiledlayout(4,3,'TileSpacing','compact','Padding','compact');

% dynamics plot, with time ranges marked
nexttile([1,3])
hold on
errorbar(dtime_avg,dtime_std)
rectangle('Position',[1 0 1 1],'LineWidth',2)
rectangle('Position',[intmin 0 (intmax-intmin) 1],'LineWidth',1)
rectangle('Position',[endmin 0 (endmax-endmin) 1],'LineWidth',1)
axis([0 500000 0 1.1])
title(strcat('Displacement Strategy Evolution, nmax=',num2str(nmax),', rep ',num2str(rep)))
xlabel('Generation')
ylabel('Probability')
hold off

% beginning plot
load('C:\Users\eschlatter\Dropbox\DispersalEvolution\files_runsims\pop_uniform_sx=2_sy=512_nbins=12.mat')
kern_beginning = mean(pop(:,1:12),1);
nexttile([2,1])
title('Beginning')
dit_beginning = kernelcosts(kern_beginning,nmax,1); %save the values for each distance

% intermediate plot
v_intermed = dtime_avg(1:intmax,:);
nexttile([2,1])
title('Middle')
weighted_avg_mid = dtimecosts_error_plot(v_intermed,nmax,intmin,1);

% ESS plot
% sample a bunch of kernels from the last 3x10^5 generations
% for each kernel sampled, calculate the components of selection
nexttile([2,1])
v_end = dtime_avg(1:endmax,:);
title('End')
weighted_avg_end = dtimecosts_error_plot(v_end,nmax,endmin,1);
legend('kernel','','mortality','kin competition','total benefit')

%
% total costs and benefits figures
%

m=16; %marker size
l=1; %errorbar line width

%beginning
nexttile([1,1])
totals_beginning = dit_beginning*kern_beginning';
hold on
plot(1,totals_beginning(1),'.','Color','#D95319','MarkerSize',m)
plot(2,totals_beginning(2),'.','Color','#EDB120','MarkerSize',m)
plot(3,totals_beginning(3),'.','Color','#7E2F8E','MarkerSize',m)
ylim([0 1])
xlim([0.5 3.5])
xticks([1 2 3])
xticklabels({'Mortality', 'Kin competition', 'Total benefit'})
hold off

%middle
nexttile([1,1])
means_mid = mean(weighted_avg_mid,1);
stds_mid = std(weighted_avg_mid);
hold on
plot(1,means_mid(1),'.','Color','#D95319','MarkerSize',m)
errorbar(1,means_mid(1),stds_mid(1),'LineWidth',l)
plot(2,means_mid(2),'.','Color','#EDB120','MarkerSize',m)
errorbar(2,means_mid(2),stds_mid(2),'Color','#EDB120','LineWidth',l)
plot(3,means_mid(3),'.','Color','#7E2F8E','MarkerSize',m)
errorbar(3,means_mid(3),stds_mid(3),'Color','#7E2F8E','LineWidth',l)
ylim([0 1])
xlim([0.5 3.5])
xticks([1 2 3])
xticklabels({'Mortality', 'Kin competition', 'Total benefit'})
hold off

%end
nexttile([1,1])
means_end = mean(weighted_avg_end,1);
stds_end = std(weighted_avg_end);
hold on
plot(1,means_end(1),'.','Color','#D95319','MarkerSize',m)
errorbar(1,means_end(1),stds_end(1),'LineWidth',l)
plot(2,means_end(2),'.','Color','#EDB120','MarkerSize',m)
errorbar(2,means_end(2),stds_end(2),'Color','#EDB120','LineWidth',l)
plot(3,means_end(3),'.','Color','#7E2F8E','MarkerSize',m)
errorbar(3,means_end(3),stds_end(3),'Color','#7E2F8E','LineWidth',l)
ylim([0 1])
xlim([0.5 3.5])
xticks([1 2 3])
xticklabels({'Mortality', 'Kin competition', 'Total benefit'})
hold off


%change dimensions of figure
fig_position = get(gcf,'Position');
fig_position(3) = 2*fig_position(3); %width
fig_position(4) = 1.5*fig_position(4); %height
set(gcf,'Position',fig_position);

%save it
exportgraphics(t,strcat('C:\Users\eschlatter\Dropbox\DispersalEvolution\output_simulations\all_20221103\comboplot2_nmax=',num2str(nmax),...
    '_rep=',num2str(rep),'.jpg'))

