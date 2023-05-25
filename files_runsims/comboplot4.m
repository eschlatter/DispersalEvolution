function comboplot4(dtime_avg,dtime_std,nmax,rep,endmin, endmax)
%plots kernel dynamics, plus costs figures for beginning and end of
%evolutionary period

figure
clf
t=tiledlayout(4,1,'Padding','compact');

% dynamics plot, with time ranges marked
ax = nexttile;
%ax.ColorOrder = (1/255)*[0 109 44; 44 162 95; 102 194 164; 178 226 226];
%ax.ColorOrder = (1/255)*[37 52 148;44 127 184;65 182 196;161 218 180; 49 163 84;0 104 55];
%ax.ColorOrder = (1/255)*[228 26 28;55 126 184;77 175 74;152 78 163;255 127 0;255 255 51];
ax.ColorOrder = (1/255)*[55 126 184;77 175 74;152 78 163;255 127 0];
ax.LineStyleOrder = {'-','--',':'};
hold on
%errorbar(dtime_avg,dtime_std,'CapSize',4)
plot(dtime_avg)
rectangle('Position',[endmin 0 (endmax-endmin) 1],'LineWidth',1)
axis([0 500000 0 1.3])
title(strcat('nmax=',num2str(nmax)', {newline}, 'Displacement Strategy Evolution'))
xlabel('Generation')
ylabel('Probability')
lgd=legend(string(0:11),'Location','North','NumColumns',6);
title(lgd,'Distance')
text(470000,1.4,'(a)','FontSize',18)
hold off

% beginning plot
load('C:\Users\eschlatter\Dropbox\DispersalEvolution\files_runsims\pop_uniform_sx=2_sy=512_nbins=12.mat')
kern_beginning = mean(pop(:,1:12),1);
nexttile()
title('Beginning Kernel')
dit_beginning = kernelcosts(kern_beginning,nmax,1); %save the values for each distance
legend('kernel','mortality','kin competition','total benefit','Location','SouthOutside','NumColumns',2)
text(12.5,1.2,'(b)','FontSize',18)

% ESS plot
% sample a bunch of kernels from the last 3x10^5 generations
% for each kernel sampled, calculate the components of selection
nexttile()
v_end = dtime_avg(1:endmax,:);
title('ESS Kernel')
weighted_avg_end = dtimecosts_error_plot(v_end,nmax,endmin,1);
text(12.5,1.2,'(c)','FontSize',18)

%
% total costs and benefits figures
%

m=24; %marker size
l=1; %errorbar line width

%beginning
nexttile()
totals_beginning = dit_beginning*kern_beginning';
means_end = mean(weighted_avg_end,1);
stds_end = std(weighted_avg_end);

hold on
plot(1,totals_beginning(1),'o','Color','#D95319','MarkerSize',m/4)
plot(2,totals_beginning(2),'o','Color','#EDB120','MarkerSize',m/4)
plot(3,totals_beginning(3),'o','Color','#7E2F8E','MarkerSize',m/4)

plot(1,means_end(1),'.','Color','#D95319','MarkerSize',m)
errorbar(1,means_end(1),stds_end(1),'Color','#D95319','LineWidth',l,'CapSize',12)
plot(2,means_end(2),'.','Color','#EDB120','MarkerSize',m)
errorbar(2,means_end(2),stds_end(2),'Color','#EDB120','LineWidth',l,'CapSize',12)
plot(3,means_end(3),'.','Color','#7E2F8E','MarkerSize',m)
errorbar(3,means_end(3),stds_end(3),'Color','#7E2F8E','LineWidth',l,'CapSize',12)

ylim([0 1])
xlim([0.5 3.5])
xticks([1 2 3])
xticklabels({'Mortality', 'Kin competition', 'Total benefit'})
ylabel('Combined cost or benefit over all distances')
title('Costs and Fitness')

qw{1} = plot(nan, 'ko','MarkerSize',m/4);
qw{2} = plot(nan, 'k.','MarkerSize',m);
legend([qw{:}], {'Beginning Kernel','ESS Kernel'}, 'location', 'north')
text(3.35,1.2,'(d)','FontSize',18)
hold off


%change dimensions of figure
fig_position = get(gcf,'Position');
fig_position(3) = 0.8*fig_position(3); %width
fig_position(4) = 3.2*fig_position(4); %height
set(gcf,'Position',fig_position);

%save it
exportgraphics(t,strcat('C:\Users\eschlatter\Dropbox\DispersalEvolution\output_simulations\all_20221103\comboplot4_nmax=',num2str(nmax),...
    '_rep=',num2str(rep),'.jpg'))

