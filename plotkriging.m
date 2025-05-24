function M=plotkriging(bsfit,bsgi,YX,X1,X2,t)

if ~exist('Figure5', 'dir')
    mkdir('Figure5')
end
for i=1:size(bsfit,1)
    data0(i)=f_num(bsfit(i,:),t);
end
for i=1:size(YX,1)
    data1(i)=f_num(YX(i,:),t);
end
data2=reshape(data1,size(X1));
figure
mesh(X1,X2,data2);
hold on;
plot3(bsgi(:,2),bsgi(:,3),data0,'.k', 'MarkerSize',10);
hold off;
saveas(gcf,['Figure5/demand_case_',num2str(t),'.jpg'])
close;
end