function comboplot(dtime_avg,dtime_std,nmax,rep, intmin, intmax, endmin, endmax)
%plots kernel dynamics, plus costs figures for beginning, middle and end of
%evolutionary period

figure
clf
t=tiledlayout(3,3);

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
nexttile([2,1])
title('Beginning')
kernelcosts(mean(pop(:,1:12),1),nmax,1);

% intermediate plot
v_intermed = dtime_avg(1:intmax,:);
nexttile([2,1])
title('Middle')
dtimecosts_error_plot(v_intermed,nmax,intmin,1);

% ESS plot
% sample a bunch of kernels from the last 3x10^5 generations
% for each kernel sampled, calculate the components of selection
nexttile([2,1])
v_end = dtime_avg(1:endmax,:);
title('End')
dtimecosts_error_plot(v_end,nmax,endmin,1);
legend('kernel','','mortality','kin competition','total benefit')

%change dimensions of figure
fig_position = get(gcf,'Position');
fig_position(3) = 2*fig_position(3);
set(gcf,'Position',fig_position);

%save it
exportgraphics(t,strcat('C:\Users\eschlatter\Dropbox\DispersalEvolution\output_simulations\all_20221103\comboplot_nmax=',num2str(nmax),...
    '_rep=',num2str(rep),'.jpg'))

