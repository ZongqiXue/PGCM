function plotcircle(x,y,r)
t= 0:pi/100:2*pi;
for i=1:size(x)
pb=patch((r*sin(t)+ x(i)),(r*cos(t)+y(i)),[0.3 0.3 1],'edgecolor','k');
%pb=patch((r*sin(t)+ x(i)),(r*cos(t)+y(i)),rand(1,3),'edgecolor','k');
alpha(pb,.3);
end
hold on
axis equal
end