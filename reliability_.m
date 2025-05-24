function f=reliability_(t,tt)
l=97441;
k=1.6395;
B=3.5;
if t<tt
    f=exp(-(t/l).^k);
else
    f=exp((B-1)*(tt/l).^k-B*(t/l).^k);
end