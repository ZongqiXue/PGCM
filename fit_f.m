function r=fit_f(data)
A=zeros(35,35);
b=zeros(35,1);
for i=1:216
    for j=1:35
        for k=1:35
        A(j,k)=A(j,k)+fit_function(i,j)*fit_function(i,k);
        end
        b(j,1)=b(j,1)+fit_function(i,j)*data(i);
    end
end
for i=270:370
    for j=1:35
        for k=1:35
        A(j,k)=A(j,k)+fit_function(i,j)*fit_function(i,k);
        end
        b(j,1)=b(j,1)+fit_function(i,j)*data(i);
    end
end
r=A\b;
end
