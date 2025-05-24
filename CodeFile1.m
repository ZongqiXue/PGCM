clc; clear all;
load Dataset1.mat
if ~exist('Figure4', 'dir')
    mkdir('Figure4')
end
  
k=1;
M=30;
for i=1:length(aname)
    figure
    eval(['a0=',aname{i},';']);
    a1=sum(a0,2);
    [a2,a3]=fit_plot_(a1,M);

    xlabel('time (hour)','FontSize',15);
    ylabel('demand (GB)','FontSize',15);
    ax=gca;
    ax.XAxis.FontSize=12;
    ax.YAxis.FontSize=12;
    %saveas(gcf,['Figure4\case',num2str(k),'.fig']);
    saveas(gcf,['Figure4\case',num2str(k),'.jpg']);
    close
    k=k+1;
end
