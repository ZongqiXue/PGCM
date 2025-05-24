clc; clear all;
load Dataset2.mat

if ~exist('Table3', 'dir')
    mkdir('Table3');
end

a0=zeros(5,6);
a1=zeros(5,6);
a2=zeros(5,6);
a3=zeros(5,6);
for case_=1:6
xlb=Case(1,case_);
xub=Case(2,case_);
yub=Case(3,case_);
ylb=Case(4,case_);   

RR1=zeros(5,200);
RR2=zeros(5,200);
term=1;
num=325;
for cof=1:0.1:1.4
%cof=1;
r=0.0039;
r=r*cof;
Ratio=foldsum(xlb,xub,ylb,yub,r,bsgi);
x=bsgi(:,2);
y=bsgi(:,3);
axis equal
[vx,vy] = voronoi(x,y);
[v,c]=voronoin([x,y]);
m=zeros(num,1);
for i=1:num
    a=c{i,1};
    for j=1:length(a)
        A_2(i,j)=0;
        A_3(i,j)=0;
        if a(j)==1
            m(i)=1;
            break
        end
    end  
    if x(i)<xlb||y(i)<ylb||x(i)>xub||y(i)>yub
        m(i)=2;
    end   
    if m(i)~=1
        A(i)=0;
        for j=1:(length(a)-1)
            A(i)=A(i)+(v(a(j),1)*v(a(j+1),2)-v(a(j+1),1)*v(a(j),2));
        end
        A(i)=A(i)+(v(a(length(a)),1)*v(a(1),2)-v(a(1),1)*v(a(length(a)),2));
        A(i)=abs(A(i)/2);
        R(i)=abs(sqrt((v(a(1),1)-x(i))^2+(v(a(1),2)-y(i))^2));
        for j=2:length(a)
            if sqrt((v(a(j),1)-x(i))^2+(v(a(j),2)-y(i))^2)>R(i)
                R(i)=abs(sqrt((v(a(j),1)-x(i))^2+(v(a(j),2)-y(i))^2));
            end
        end
       for j=1:(length(a)-1)
          for k=1:num
              if k~=i
                  b=c{k,1};
                  for p=1:length(b)
                      if p==1
                          if b(p)==a(j)  
                              if b(length(b))==a(j+1)||b(p+1)==a(j+1)
                                  n(i,j)=k;
                                  break
                              end
                          end
                      elseif p==length(b)                          
                          if b(p)==a(j)                          
                              if b(p-1)==a(j+1)||b(1)==a(j+1)                            
                                  n(i,j)=k;                              
                                  break                         
                              end                             
                          end
                      else
                          if b(p)==a(j)
                              if b(p-1)==a(j+1)||b(p+1)==a(j+1)
                                  n(i,j)=k;
                                  break
                              end
                          end
                      end
                  end
              end
          end
       end
       for k=1:num
              if k~=i
                  b=c{k,1};
                  for p=1:length(b)
                      if p==1
                          if b(1)==a(length(a))  
                              if b(length(b))==a(1)||b(p+1)==a(1)
                                  n(i,length(a))=k;
                                  break
                              end
                          end
                      elseif p==length(b)                          
                          if b(p)==a(length(a))                          
                              if b(p-1)==a(1)||b(1)==a(1)                            
                                  n(i,length(a))=k;                              
                                  break                         
                              end                             
                          end
                      else
                          if b(p)==a(length(a))
                              if b(p-1)==a(1)||b(p+1)==a(1)
                                  n(i,length(a))=k;
                                  break
                              end
                          end
                      end
                  end
              end
        end
    end
end

for i=1:num
    if m(i)==0
       a=c{i,1};
       for j=1:length(a)
           if sqrt((x(i)-x(n(i,j)))^2+(y(i)-y(n(i,j)))^2)<r
           k=-(x(i)-x(n(i,j)))/(y(i)-y(n(i,j)));
           b=y(n(i,j))-(y(i)+y(n(i,j)))/2+k*(x(i)+x(n(i,j)))/2;
           ro=roots([(1+k^2) (-2*x(n(i,j))-2*k*b) (b^2+x(n(i,j))^2-r^2)]);
           x_1=ro(1,1);
           y_1=k*x_1+(y(i)+y(n(i,j)))/2-k*(x(i)+x(n(i,j)))/2;
           x_2=ro(2,1);
           y_2=k*x_2+(y(i)+y(n(i,j)))/2-k*(x(i)+x(n(i,j)))/2;
           r1=[x_1-x(i),y_1-y(i)];
           r2=[x_2-x(i),y_2-y(i)];
           alpha=acos(dot(r1,r2)/(norm(r1)*norm(r2)));
           S_tri=abs(x_1*y_2-x_2*y_1+x_2*y(i)-x(i)*y_2+x(i)*y_1-x_1*y(i))/2;
           A_2(i,j)=alpha*r^2-S_tri*2;
           xx_2(i,j)=(x(i)+x(n(i,j)))/2;
           yy_2(i,j)=(y(i)+y(n(i,j)))/2;
           end
       end
    end
end
for i=1:num
    if m(i)==0
        a=c{i,1};
        for j=1:length(a)
            if sqrt((x(i)-x(n(i,j)))^2+(y(i)-y(n(i,j)))^2)<r
            if j==length(a)
                if sqrt((x(i)-x(n(i,1)))^2+(y(i)-y(n(i,1)))^2)<r&sqrt((x(n(i,1))-x(n(i,j)))^2+(y(n(i,1))-y(n(i,j)))^2)<r
                xx_3(i,j)=(x(i)+x(n(i,j))+x(n(i,1)))/3;
                yy_3(i,j)=(y(i)+y(n(i,j))+y(n(i,1)))/3;
                sign1=1;
                k=-(x(i)-x(n(i,j)))/(y(i)-y(n(i,j)));
                b=y(n(i,j))-(y(i)+y(n(i,j)))/2+k*(x(i)+x(n(i,j)))/2;
                ro=roots([(1+k^2) (-2*x(n(i,j))-2*k*b) (b^2+x(n(i,j))^2-r^2)]);
                sign11=(ro(1,1)-x(n(i,1)))^2+(k*ro(1,1)+(y(i)+y(n(i,j)))/2-k*(x(i)+x(n(i,j)))/2-y(n(i,1)))^2<=r^2;
                sign12=(ro(2,1)-x(n(i,1)))^2+(k*ro(2,1)+(y(i)+y(n(i,j)))/2-k*(x(i)+x(n(i,j)))/2-y(n(i,1)))^2<=r^2;
                if sign11&&sign12
                    sign1=0;
                elseif sign11==1&&sign12==0
                    x_1=ro(1,1);
                    y_1=k*ro(1,1)+(y(i)+y(n(i,j)))/2-k*(x(i)+x(n(i,j)))/2;
                elseif sign11==0&&sign12==1
                    x_1=ro(2,1);
                    y_1=k*ro(2,1)+(y(i)+y(n(i,j)))/2-k*(x(i)+x(n(i,j)))/2;
                end
                sign2=1;
                k=-(x(i)-x(n(i,1)))/(y(i)-y(n(i,1)));
                b=y(n(i,1))-(y(i)+y(n(i,1)))/2+k*(x(i)+x(n(i,1)))/2;
                ro=roots([(1+k^2) (-2*x(n(i,1))-2*k*b) (b^2+x(n(i,1))^2-r^2)]);
                sign21=(ro(1,1)-x(n(i,j)))^2+(k*ro(1,1)+(y(i)+y(n(i,1)))/2-k*(x(i)+x(n(i,1)))/2-y(n(i,j)))^2<=r^2;
                sign22=(ro(2,1)-x(n(i,j)))^2+(k*ro(2,1)+(y(i)+y(n(i,1)))/2-k*(x(i)+x(n(i,1)))/2-y(n(i,j)))^2<=r^2;
                if sign21&&sign22
                    sign2=0;
                elseif sign21==1&&sign22==0
                    x_2=ro(1,1);
                    y_2=k*ro(1,1)+(y(i)+y(n(i,1)))/2-k*(x(i)+x(n(i,1)))/2;
                elseif sign21==0&&sign22==1
                    x_2=ro(2,1);
                    y_2=k*ro(2,1)+(y(i)+y(n(i,1)))/2-k*(x(i)+x(n(i,1)))/2;
                end
                sign3=1;
                k=-(x(n(i,j))-x(n(i,1)))/(y(n(i,j))-y(n(i,1)));
                b=y(n(i,j))-(y(n(i,1))+y(n(i,j)))/2+k*(x(n(i,1))+x(n(i,j)))/2;
                ro=roots([(1+k^2) (-2*x(n(i,j))-2*k*b) (b^2+x(n(i,j))^2-r^2)]);
                sign31=(ro(1,1)-x(i))^2+(k*ro(1,1)+(y(n(i,1))+y(n(i,j)))/2-k*(x(n(i,1))+x(n(i,j)))/2-y(i))^2<=r^2;
                sign32=(ro(2,1)-x(i))^2+(k*ro(2,1)+(y(n(i,1))+y(n(i,j)))/2-k*(x(n(i,1))+x(n(i,j)))/2-y(i))^2<=r^2;
                if sign31&&sign32
                    sign3=0;
                    m1=ro(1,1);
                    m2=ro(2,1);
                    n1=k*ro(1,1)+(y(n(i,1))+y(n(i,j)))/2-k*(x(n(i,1))+x(n(i,j)))/2;
                    n2=k*ro(2,1)+(y(n(i,1))+y(n(i,j)))/2-k*(x(n(i,1))+x(n(i,j)))/2;
                elseif sign31==1&&sign32==0
                    x_3=ro(1,1);
                    y_3=k*ro(1,1)+(y(n(i,1))+y(n(i,j)))/2-k*(x(n(i,1))+x(n(i,j)))/2;
                elseif sign31==0&&sign32==1
                    x_3=ro(2,1);
                    y_3=k*ro(2,1)+(y(n(i,1))+y(n(i,j)))/2-k*(x(n(i,1))+x(n(i,j)))/2;
                end
                if sign1==0
                    A_3(i,j)=A_2(i,j);
                elseif sign2==0
                    A_3(i,j)=A_2(i,j+1);
                elseif sign3==0                              
                    r1=[m1-x(n(i,j)),n1-y(n(i,j))];       
                    r2=[m2-x(n(i,j)),n2-y(n(i,j))];                     
                    alpha=acos(dot(r1,r2)/(norm(r1)*norm(r2)));                              
                    S_tri=abs(m1*n2-m2*n1+m2*y(n(i,j))-x(n(i,j))*n2+x(n(i,j))*n1-m1*y(n(i,j)))/2;
                    A_3(i,j)=alpha*r^2-S_tri*2;
                else                  
                    r1=[x_1-x(i),y_1-y(i)];         
                    r2=[x_2-x(i),y_2-y(i)];
                    alpha=acos(dot(r1,r2)/(norm(r1)*norm(r2)));                              
                    s1=alpha*r^2/2-abs(x_1*y_2-x_2*y_1+x_2*y(i)-x(i)*y_2+x(i)*y_1-x_1*y(i))/2;
                    r1=[x_1-x(n(i,j)),y_1-y(n(i,j))];         
                    r2=[x_3-x(n(i,j)),y_3-y(n(i,j))];
                    alpha=acos(dot(r1,r2)/(norm(r1)*norm(r2)));                              
                    s2=alpha*r^2/2-abs(x_1*y_3-x_3*y_1+x_3*y(n(i,j))-x(n(i,j))*y_3+x(n(i,j))*y_1-x_1*y(n(i,j)))/2;
                    r1=[x_2-x(n(i,1)),y_2-y(n(i,1))];         
                    r2=[x_3-x(n(i,1)),y_3-y(n(i,1))];
                    alpha=acos(dot(r1,r2)/(norm(r1)*norm(r2)));                              
                    s3=alpha*r^2/2-abs(x_3*y_2-x_2*y_3+x_2*y(n(i,1))-x(n(i,1))*y_2+x(n(i,1))*y_3-x_3*y(n(i,1)))/2;
                    A_3(i,j)=s1+s2+s3+abs(x_1*y_2-x_2*y_1+x_2*y_3-x_3*y_2+x_3*y_1-x_1*y_3)/2;
                    %sign_(i,j)=1;
                    %x1_(i,j)=x_1;
                    %y1_(i,j)=y_1;
                    %x2_(i,j)=x_2;
                    %y2_(i,j)=y_2;
                    %x3_(i,j)=x_3;
                    %y3_(i,j)=y_3;
                end 
                end
            else
                if sqrt((x(i)-x(n(i,j+1)))^2+(y(i)-y(n(i,j+1)))^2)<r&sqrt((x(n(i,j+1))-x(n(i,j)))^2+(y(n(i,j+1))-y(n(i,j)))^2)<r
                xx_3(i,j)=(x(i)+x(n(i,j))+x(n(i,j+1)))/3;
                yy_3(i,j)=(y(i)+y(n(i,j))+y(n(i,j+1)))/3;
                sign1=1;
                k=-(x(i)-x(n(i,j)))/(y(i)-y(n(i,j)));
                b=y(n(i,j))-(y(i)+y(n(i,j)))/2+k*(x(i)+x(n(i,j)))/2;
                ro=roots([(1+k^2) (-2*x(n(i,j))-2*k*b) (b^2+x(n(i,j))^2-r^2)]);
                sign11=(ro(1,1)-x(n(i,j+1)))^2+(k*ro(1,1)+(y(i)+y(n(i,j)))/2-k*(x(i)+x(n(i,j)))/2-y(n(i,j+1)))^2<=r^2;
                sign12=(ro(2,1)-x(n(i,j+1)))^2+(k*ro(2,1)+(y(i)+y(n(i,j)))/2-k*(x(i)+x(n(i,j)))/2-y(n(i,j+1)))^2<=r^2;
                if sign11&&sign12
                    sign1=0;
                elseif sign11==1&&sign12==0
                    x_1=ro(1,1);
                    y_1=k*ro(1,1)+(y(i)+y(n(i,j)))/2-k*(x(i)+x(n(i,j)))/2;
                elseif sign11==0&&sign12==1
                    x_1=ro(2,1);
                    y_1=k*ro(2,1)+(y(i)+y(n(i,j)))/2-k*(x(i)+x(n(i,j)))/2;
                end
                sign2=1;
                k=-(x(i)-x(n(i,j+1)))/(y(i)-y(n(i,j+1)));
                b=y(n(i,j+1))-(y(i)+y(n(i,j+1)))/2+k*(x(i)+x(n(i,j+1)))/2;
                ro=roots([(1+k^2) (-2*x(n(i,j+1))-2*k*b) (b^2+x(n(i,j+1))^2-r^2)]);
                sign21=(ro(1,1)-x(n(i,j)))^2+(k*ro(1,1)+(y(i)+y(n(i,j+1)))/2-k*(x(i)+x(n(i,j+1)))/2-y(n(i,j)))^2<=r^2;
                sign22=(ro(2,1)-x(n(i,j)))^2+(k*ro(2,1)+(y(i)+y(n(i,j+1)))/2-k*(x(i)+x(n(i,j+1)))/2-y(n(i,j)))^2<=r^2;
                if sign21&&sign22
                    sign2=0;
                elseif sign21==1&&sign22==0
                    x_2=ro(1,1);
                    y_2=k*ro(1,1)+(y(i)+y(n(i,j+1)))/2-k*(x(i)+x(n(i,j+1)))/2;
                elseif sign21==0&&sign22==1
                    x_2=ro(2,1);
                    y_2=k*ro(2,1)+(y(i)+y(n(i,j+1)))/2-k*(x(i)+x(n(i,j+1)))/2;
                end
                sign3=1;
                k=-(x(n(i,j))-x(n(i,j+1)))/(y(n(i,j))-y(n(i,j+1)));
                b=y(n(i,j))-(y(n(i,j+1))+y(n(i,j)))/2+k*(x(n(i,j+1))+x(n(i,j)))/2;
                ro=roots([(1+k^2) (-2*x(n(i,j))-2*k*b) (b^2+x(n(i,j))^2-r^2)]);
                sign31=(ro(1,1)-x(i))^2+(k*ro(1,1)+(y(n(i,j+1))+y(n(i,j)))/2-k*(x(n(i,j+1))+x(n(i,j)))/2-y(i))^2<=r^2;
                sign32=(ro(2,1)-x(i))^2+(k*ro(2,1)+(y(n(i,j+1))+y(n(i,j)))/2-k*(x(n(i,j+1))+x(n(i,j)))/2-y(i))^2<=r^2;
                if sign31&&sign32
                    sign3=0;
                    m1=ro(1,1);
                    m2=ro(2,1);
                    n1=k*ro(1,1)+(y(n(i,j+1))+y(n(i,j)))/2-k*(x(n(i,j+1))+x(n(i,j)))/2;
                    n2=k*ro(2,1)+(y(n(i,j+1))+y(n(i,j)))/2-k*(x(n(i,j+1))+x(n(i,j)))/2;
                elseif sign31==1&&sign32==0
                    x_3=ro(1,1);
                    y_3=k*ro(1,1)+(y(n(i,j+1))+y(n(i,j)))/2-k*(x(n(i,j+1))+x(n(i,j)))/2;
                elseif sign31==0&&sign32==1
                    x_3=ro(2,1);
                    y_3=k*ro(2,1)+(y(n(i,j+1))+y(n(i,j)))/2-k*(x(n(i,j+1))+x(n(i,j)))/2;
                end
                if sign1==0
                    A_3(i,j)=A_2(i,j);
                elseif sign2==0
                    A_3(i,j)=A_2(i,j+1);
                elseif sign3==0                              
                    r1=[m1-x(n(i,j)),n1-y(n(i,j))];       
                    r2=[m2-x(n(i,j)),n2-y(n(i,j))];                     
                    alpha=acos(dot(r1,r2)/(norm(r1)*norm(r2)));                              
                    S_tri=abs(m1*n2-m2*n1+m2*y(n(i,j))-x(n(i,j))*n2+x(n(i,j))*n1-m1*y(n(i,j)))/2;
                    A_3(i,j)=alpha*r^2-S_tri*2;
                else                  
                    r1=[x_1-x(i),y_1-y(i)];         
                    r2=[x_2-x(i),y_2-y(i)];
                    alpha=acos(dot(r1,r2)/(norm(r1)*norm(r2)));                              
                    s1=alpha*r^2/2-abs(x_1*y_2-x_2*y_1+x_2*y(i)-x(i)*y_2+x(i)*y_1-x_1*y(i))/2;
                    %s1_(i,j)=s1
                    %alpha1(i,j)=alpha;
                    %result(i,j)=dot(r1,r2)/(norm(r1)*norm(r2));
                    r1=[x_1-x(n(i,j)),y_1-y(n(i,j))];         
                    r2=[x_3-x(n(i,j)),y_3-y(n(i,j))];
                    alpha=acos(dot(r1,r2)/(norm(r1)*norm(r2)));                              
                    s2=alpha*r^2/2-abs(x_1*y_3-x_3*y_1+x_3*y(n(i,j))-x(n(i,j))*y_3+x(n(i,j))*y_1-x_1*y(n(i,j)))/2;
                    %alpha2(i,j)=alpha;
                    r1=[x_2-x(n(i,j+1)),y_2-y(n(i,j+1))];         
                    r2=[x_3-x(n(i,j+1)),y_3-y(n(i,j+1))];
                    alpha=acos(dot(r1,r2)/(norm(r1)*norm(r2)));                              
                    s3=alpha*r^2/2-abs(x_3*y_2-x_2*y_3+x_2*y(n(i,j+1))-x(n(i,j+1))*y_2+x(n(i,j+1))*y_3-x_3*y(n(i,j+1)))/2;                 
                    %alpha3(i,j)=alpha
                    A_3(i,j)=s1+s2+s3+abs(x_1*y_2-x_2*y_1+x_2*y_3-x_3*y_2+x_3*y_1-x_1*y_3)/2;
                    %sign_(i,j)=1;
                    %x1_(i,j)=x_1;
                    %y1_(i,j)=y_1;
                    %x2_(i,j)=x_2;
                    %y2_(i,j)=y_2;
                    %x3_(i,j)=x_3;
                    %y3_(i,j)=y_3;
                end
                end
            end
            end
        end
    end
end
for i=1:num
    if m(i)==0
        A_1(i)=pi*r^2;    
        a=c{i,1};
        for j=1:length(a)
            A_1(i)=A_1(i)-A_2(i,j)+A_3(i,j);
        end
    end
end
clear R1 R2
a3(round(10*(cof-1)+1),case_)=sum(sum(A_3))/3;
a2(round(10*(cof-1)+1),case_)=(sum(sum(A_2))-sum(sum(A_3)))/2;
a1(round(10*(cof-1)+1),case_)=length(find(m==0))*pi*r^2-sum(sum(A_2));
a0(round(10*(cof-1)+1),case_)=length(find(m==0))*pi*r^2;
tt=0;
for t=1:50:10000
num=length(A);
tt=tt+1;
R1(tt)=0;
k=0;
for i=1:num
    if m(i)==0
       k=k+1;
       R_(i)=0;
       a=c{i,1};
       for j=1:length(a)
           R_(i)=R_(i)+abs(A_2(i,j)*reliability(t)*f_num(zz_2(i,j,:),t)-A_3(i,j)*reliability(t)*reliability(t)*f_num(zz_3(i,j,:),t));
       end
       R_(i)=R_(i)/(R_(i)+max([A_1(i),0])*f_num(bsfit(i,:),t));
       if R_(i)<inf
       R1(tt)=R1(tt)+R_(i);
       end
    end 
end
R1(tt)=R1(tt)/k;
end
RR1(term,:)=R1;
tt=0;
for t=1:50:10000
num=length(A);
tt=tt+1;
R2(tt)=0;
k=0;
R_=0;
R__=0;
for i=1:num
    if m(i)==0
       k=k+1;
       a=c{i,1};
       rr=0;
       rrr=0;
       for j=1:length(a)
           rr=rr+(-A_2(i,j)*(1-reliability(t)^2)*f_num(zz_2(i,j,:),t)/2+A_3(i,j)*(1-reliability(t)^3)*f_num(zz_3(i,j,:),t)/3);
           rrr=rrr+(-A_2(i,j)*f_num(zz_2(i,j,:),t)/2+A_3(i,j)*f_num(zz_3(i,j,:),t)/3);
       end
       R_=R_+abs(pi*r*r*(1-reliability(t))*f_num(bsfit(i,:),t)+rr);
       R__=R__+abs(pi*r*r*f_num(bsfit(i,:),t)+rrr);
    end    
end
R2(tt)=1-R_/R__;
end

RR2(term,:)=(1-Ratio(1))*R2;
term=term+1;
end
r=0.0045;
angle=0:pi/100:2*pi;
figure
for i=1:num
    pb=patch((r*sin(angle)+ x(i)),(r*cos(angle)+y(i)),[0.4*rand(1,1) 1-0.5*rand(1,1) 1-0.5*rand(1,1)],'facealpha',0.3,'edgecolor','k');    
    hold on 
    scatter(x(i),y(i),'.','k');
end
axis([xlb,xub,ylb,yub])
axis equal
xlim([xlb,xub]);
ylim([ylb,yub]);
xlabel('longitude ');
ylabel('latitude ');
box on;
set(gca, 'FontSize', 15);
saveas(gcf,['Table3\point[',num2str(xlb),',',num2str(xub),',',num2str(ylb),',',num2str(yub),'].jpg']);
close
t=1:50:10000;
figure
hold on
for i=1:5
    plot(t,RR1(i,:));
end
legend('r=330m','r=365m','r=400m','r=435m','r=470m','NumColumns',2);
xlabel('t/hour');
ylabel('resilience');
%title('resilience');
ylim([min(min(RR1))-0.035,max(max(RR1))+0.5*(max(max(RR1))-min(min(RR1)))]);
yticks(floor(10*(min(min(RR1))-0.035))/10:0.1:0.1+min(1,floor(10*max(max(RR1))+0.035)/10));
box on;
set(gca, 'FontSize', 15);
saveas(gcf,['Table3\resilience[',num2str(xlb),',',num2str(xub),',',num2str(ylb),',',num2str(yub),'].jpg']);
close
figure
hold on
for i=1:5
    plot(t,RR2(i,:));
end
legend('r=330m','r=365m','r=400m','r=435m','r=470m','NumColumns',2);
xlabel('t/hour');
ylabel('reliability');
%title('reliability');
ylim([min(min(RR2))-0.001,max(max(RR2))+0.35*(max(max(RR2))-min(min(RR2)))]);
%yticks(floor(10*(min(min(RR2))-0.001))/10:0.1:1);
box on;
set(gca, 'FontSize', 15);
saveas(gcf,['Table3\reliability[',num2str(xlb),',',num2str(xub),',',num2str(ylb),',',num2str(yub),'].jpg']);
close
end


