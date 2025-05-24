function [x]=variogram(h,m,n)
c=m;
a=0.03;
if n==1
    if h<a
        x=c*(1.5*h/a-0.5*(h/a)^3);
    else
        x=c;
    end
elseif n==2
    x=c*(1-exp(-3*h/a));
else
    x=c*(1-exp(-3*(h/a)^2));
end
end

        