function Ratio=foldsum(xlb,xub,ylb,yub,r,bsgi)
x=bsgi(:,2);
y=bsgi(:,3);
num=size(x,1);
angle=linspace(0, 2*pi, 100);
filename='testpicture.jpg';
figure
hold on
for i=1:num
    pb=patch((r*sin(angle)+ x(i)),(r*cos(angle)+y(i)),[0 1 0],'facealpha',0.2,'edgecolor','none');    
    hold on 
end
axis off
axis equal
xlim([xlb,xub]);
ylim([ylb,yub]);
saveas(gcf,filename);
close

img=imread(filename);
img=img(50:585,185:720,3);

Ratio=zeros(5,1);
step=[220,177,138,103];
Ratio(1)=length(find(img>=step(1)));
for i=2:4
    Ratio(i)=length(find(img>=step(i) & img<step(i-1)));
end
Ratio(5)=length(find(img<step(4)));
Ratio=Ratio/sum(Ratio);
delete(filename);
end