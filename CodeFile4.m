clc; clear all;
load Dataset3.mat
if ~exist('Table4', 'dir')
    mkdir('Table4');
end
for case_=1:6

xlb=Case(1,case_);
xub=Case(2,case_);
yub=Case(3,case_);
ylb=Case(4,case_);   
Ratio=foldsum(xlb,xub,ylb,yub,r,bsgi);
RRR1=reshape(Resilience(case_,:,:),4,200);
RRR2=(1-Ratio(1))*reshape(Reliability(case_,:,:),4,200);

figure
t=1:50:10000;
hold on
plot(t,RRR1(3,:),'c');
plot(t,RRR1(4,:),'c');
region=fill([t fliplr(t)],[RRR1(3,:) fliplr(RRR1(4,:))],'c','FaceAlpha',0.25,'EdgeColor','none');
line2=plot(t,RRR1(2,:),'b');
line1=plot(t,RRR1(1,:),'r');
legend([line1 line2 region],'resilience without alarm','mean resilience with alarm','95% confidence interval');
xlabel('t/hour');
ylabel('resilience');
%title('resilience');
ylim([min(min(RRR1))-0.0005,max(max(RRR1))+0.32*(max(max(RRR1))-min(min(RRR1)))+0.003]);
box on;
set(gca, 'FontSize', 15);
saveas(gcf,['Table4\resilience_',num2str(xlb),',',num2str(xub),',',num2str(ylb),',',num2str(yub),'].jpg']);
close

figure
t=1:50:10000;
hold on
plot(t,RRR2(3,:),'c');
plot(t,RRR2(4,:),'c');
region=fill([t fliplr(t)],[RRR2(3,:) fliplr(RRR2(4,:))],'c','FaceAlpha',0.3,'EdgeColor','none');
line2=plot(t,RRR2(2,:),'b');
line1=plot(t,RRR2(1,:),'r');
legend([line1 line2 region],'reliability without alarm','mean reliability with alarm','95% confidence interval');
xlabel('t/hour');
ylabel('reliability');
%title('reliability');
ylim([min(min(RRR2))-0.015,max(max(RRR2))+0.15*(max(max(RRR2))-min(min(RRR2)))+0.025]);
box on;
set(gca, 'FontSize', 15);
saveas(gcf,['Table4\reliability_',num2str(tini),'[',num2str(xlb),',',num2str(xub),',',num2str(ylb),',',num2str(yub),'].jpg']);
close
end