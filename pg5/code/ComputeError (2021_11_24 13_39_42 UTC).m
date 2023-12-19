function [L1,Linf,img]=ComputeError(U,element,node,Coeff,C,u)
[N,~]=size(element);
x=[-sqrt(3/5),0,sqrt(3/5)];
A=[5/9,8/9,5/9];
B=[0,0,1;
    0,1,1;
    1,0,1];
L1=zeros(1,N);
Linf=zeros(1,N);
img=zeros(9*N,3);
for m=1:N
    P=element(m,2:4);
    Q=B\node(P+1,2:3);
    for s=1:3
        for k=1:3
            xx=Q(1,1)*(1+x(k))/2+Q(2,1)*(1-x(k))*(1+x(s))/4+Q(3,1);
            yy=Q(1,2)*(1+x(k))/2+Q(2,2)*(1-x(k))*(1+x(s))/4+Q(3,2);
            tep=u(xx,yy)-Compute_Uh_2D(xx,yy,U,element,node,Coeff,C);
            L1(m)=L1(m)+A(s)*A(k)*abs(tep)*(1-x(k))/8;
            img(9*(m-1)+3*(s-1)+k,1)=xx;
            img(9*(m-1)+3*(s-1)+k,2)=yy;
            img(9*(m-1)+3*(s-1)+k,3)=tep;
        end
    end
    L1(m)=L1(m)*abs(Q(1,1)*Q(2,2)-Q(2,1)*Q(1,2));
end
L1=sum(L1);
Linf=max(abs(img(:,3)));
end