function U=ComputeCoef(N,d,c,f)
h=1/N;
xj=zeros(N+1,1);
for j=1:N+1
    xj(j)=(j-1)*h;
end
Index=zeros(N,2);
for e=2:N-1
    Index(e,1)=e-1;
    Index(e,2)=e;
end
Index(1,1)=0;
Index(1,2)=1;
Index(N,1)=N-1;
Index(N,2)=0;
%Assemble
K=zeros(N-1,N-1);
for e=2:N-1
    K(Index(e,1),Index(e,1))=K(Index(e,1),Index(e,1))+Gauss3(@(x) d(x)+c(x)*(xj(e+1)-x)^2,xj(e),xj(e+1));
    K(Index(e,2),Index(e,2))=K(Index(e,2),Index(e,2))+Gauss3(@(x) d(x)+c(x)*(x-xj(e))^2,xj(e),xj(e+1));
    K(Index(e,1),Index(e,2))=K(Index(e,1),Index(e,2))+Gauss3(@(x) -d(x)+c(x)*(xj(e+1)-x)*(x-xj(e)),xj(e),xj(e+1));
    K(Index(e,2),Index(e,1))=K(Index(e,2),Index(e,1))+Gauss3(@(x) -d(x)+c(x)*(xj(e+1)-x)*(x-xj(e)),xj(e),xj(e+1));
end
K(1,1)=K(1,1)+Gauss3(@(x) d(x)+c(x)*(x-xj(1))^2,xj(1),xj(2));
K(N-1,N-1)=K(N-1,N-1)+Gauss3(@(x) d(x)+c(x)*(xj(N+1)-x)^2,xj(N),xj(N+1));
K=K./h^2;
%Calculate F
F=zeros(N-1,1);
for e=2:N-1
    F(Index(e,1))=F(Index(e,1))+Gauss3(@(x) f(x)*(xj(e+1)-x)/h,xj(e),xj(e+1));
    F(Index(e,2))=F(Index(e,2))+Gauss3(@(x) f(x)*(x-xj(e))/h,xj(e),xj(e+1));
end
F(1)=F(1)+Gauss3(@(x) f(x)*(x-xj(1))/h,xj(1),xj(2));
F(N-1)=F(N-1)+Gauss3(@(x) f(x)*(xj(N+1)-x)/h,xj(N),xj(N+1));
%solve the equation
U=K\F;
end