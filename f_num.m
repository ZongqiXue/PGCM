function m=f_num(r,t)
m=0;
for i=1:35
    m=m+r(i)*fit_function(t,i);
end
if m<0
    m=0;
end
end