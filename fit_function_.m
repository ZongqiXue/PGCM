function a=fit_function_(t,c,M)

if c==1
    a=1;
elseif c<=M+1
    a=sin(2*pi*t*(c-1)/168);
else
    a=cos(2*pi*t*(c-M-1)/168);
end
end
