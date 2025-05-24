function [f,g]=fit_plot_(data,M)
warning('off');
t=1:370;
n=fit_f(data);
%n=fit_f_(data,M);
f=zeros(370,1);
for i=1:370
    %f(i)=f_num_(n,i,M);
    f(i)=f_num(n,i);
end
n=fit_f_(data,3*M);
g=zeros(370,1);
for i=1:370
    g(i)=f_num_(n,i,3*M);
end
warning('on');
plot(t,f,'c-');
hold on
plot(t,g,'g-');
plot(t,data(1:370),'r-');
%legend('fitted demand','original demand data');
%saveas(gcf,[m,'.jpg']);
end

    