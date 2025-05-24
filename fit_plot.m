function f=fit_plot(data)
t=1:370;
n=fit_f(data);
f=zeros(370,1);
for i=1:370
    f(i)=f_num(n,i);
end
plot(t,f_num(n,t),'c-');
hold on
plot(t,data(1:370),'r-');
%legend('fitted demand','original demand data');
%saveas(gcf,[m,'.jpg']);
end

    