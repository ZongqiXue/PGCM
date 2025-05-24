function r=fit_f_(data,M)
A=zeros(2*M+1,2*M+1);
b=zeros(2*M+1,1);
for i=1:370
    for j=1:2*M+1
        for k=1:2*M+1
        A(j,k)=A(j,k)+fit_function_(i,j,M)*fit_function_(i,k,M);
        end
        b(j,1)=b(j,1)+fit_function_(i,j,M)*data(i);
    end
end
r=A\b;
end
