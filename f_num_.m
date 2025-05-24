function m=f_num_(r,t,M)
m=0;
for i=1:2*M+1
    m=m+r(i)*fit_function_(t,i,M);
end
if m<0
    m=0;
end
end