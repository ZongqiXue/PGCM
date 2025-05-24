function [r R0 R1 R2]=holistic_measure(bsgi,bsfit,xx_2,yy_2,zz_2,xx_3,yy_3,zz_3,m,n,c,xub,xlb,yub,ylb)
num=325;
x=bsgi(:,2);
y=bsgi(:,3);
axis equal
[vx,vy] = voronoi(x,y);
[v,c]=voronoin([x,y]);
for i=1:num
    a=c{i,1};
    m(i)=0;
    for j=1:length(a)
        if a(j)==1
            m(i)=1;
            break
        end
    end
    if m(i)==0
        for j=1:length(a)
            if v(a(j),1)<xlb||v(a(j),2)<ylb||v(a(j),1)>xub||v(a(j),2)>yub
                m(i)=2;
                break
            end
        end
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
r=0;
r_=0;
for i=1:num
    if m(i)~=1
        if r_<R(i)
            r_=R(i);
        end
    end
    if m(i)==0
        if r<R(i)
            r=R(i);
        end
    end
end
for i=1:num
    if m(i)==0
       a=c{i,1};
       for j=1:length(a)
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
for i=1:num
    if m(i)==0
        a=c{i,1};
        for j=1:length(a)
            if j==length(a)
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
            else
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
for i=1:num
    if m(i)==0
        A_1(i)=pi*r^2;    
        a=c{i,1};
        for j=1:length(a)
            A_1(i)=A_1(i)-A_2(i,j)+A_3(i,j);
        end
    end
end
for t=1:170
num=length(A);
R0(t)=0;
k=0;
for i=1:num
    if m(i)==0
       k=k+1;
       R_(i)=0;
       a=c{i,1};
       for j=1:length(a)
           R_(i)=R_(i)+abs(A_2(i,j)*f_num(zz_2(i,j,:),t)-A_3(i,j)*f_num(zz_3(i,j,:),t));
       end
       R_(i)=R_(i)/(R_(i)+abs(A_1(i))*f_num(bsfit(i,:),t));
       R0(t)=R0(t)+R_(i);
    end 
end
R0(t)=R0(t)/k;
end
tt=0;
for t=1:50:35000
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
       R_(i)=R_(i)/(R_(i)+abs(A_1(i))*f_num(bsfit(i,:),t));
       R1(tt)=R1(tt)+R_(i);
    end 
end
R1(tt)=R1(tt)/k;
end
tt=0;
for t=1:50:3000
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
       for j=1:length(a)
           rr=rr+(-A_2(i,j)*reliability(t)*f_num(zz_2(i,j,:),t)+A_3(i,j)*reliability(t)*reliability(t)*f_num(zz_3(i,j,:),t));
       end
       R_=R_+abs(pi*rr*rr*f_num(bsfit(i,:),t)+rr);
       R__=R__+A_1(i)*f_num(bsfit(i,:),t)+A_2(i,j)*f_num(zz_2(i,j,:),t)/2+A_3(i,j)*f_num(zz_3(i,j,:),t)/3;
    end    
end
R2(tt)=1-(1-reliability(t))*R_/R__;
end
figure
plot(R0)
xlabel('t/hour');
ylabel('resilience');
title('resilience');
figure
t=1:50:35000;
plot(t,R1)
xlabel('t/hour');
ylabel('resilience');
title('resilience');
figure
t=1:50:3000;
plot(t,R2)
xlabel('t/hour');
ylabel('reliability');
title('reliability');
end