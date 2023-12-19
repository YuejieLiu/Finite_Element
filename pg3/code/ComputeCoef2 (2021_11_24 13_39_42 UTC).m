function U=ComputeCoef2(N,d,c,f)
h=1/N;
xj=zeros(N+1,1);
for j=1:N+1
    xj(j)=(j-1)*h;
end
Index=zeros(N,3);
for e=2:N-1
    Index(e,1)=e-1;
    Index(e,2)=e;
    Index(e,3)=N+e-1;
end
Index(1,1)=0;
Index(1,2)=1;
Index(1,3)=N;
Index(N,1)=N-1;
Index(N,2)=0;
Index(N,3)=2*N-1;
%Assemble
K=zeros(2*N-1,2*N-1);
for e=2:N-1
    K(Index(e,1),Index(e,1))=K(Index(e,1),Index(e,1))+Gauss3(@(x) d(x)*(4*x-3*xj(e+1)-xj(e))^2+c(x)*((2*x-xj(e)-xj(e+1))*(x-xj(e+1)))^2,xj(e),xj(e+1));
    K(Index(e,2),Index(e,2))=K(Index(e,2),Index(e,2))+Gauss3(@(x) d(x)*(4*x-3*xj(e)-xj(e+1))^2+c(x)*((2*x-xj(e+1)-xj(e))*(x-xj(e)))^2,xj(e),xj(e+1));
    K(Index(e,3),Index(e,3))=K(Index(e,3),Index(e,3))+Gauss3(@(x) d(x)*(4*(xj(e+1)-2*x+xj(e)))^2+c(x)*(4*(x-xj(e))*(xj(e+1)-x))^2,xj(e),xj(e+1));
    K(Index(e,1),Index(e,2))=K(Index(e,1),Index(e,2))+Gauss3(@(x) d(x)*(4*x-3*xj(e+1)-xj(e))*(4*x-3*xj(e)-xj(e+1))+...
        c(x)*((2*x-xj(e)-xj(e+1))*(x-xj(e+1)))*((2*x-xj(e+1)-xj(e))*(x-xj(e))),xj(e),xj(e+1));
    K(Index(e,2),Index(e,1))=K(Index(e,2),Index(e,1))+Gauss3(@(x) d(x)*(4*x-3*xj(e+1)-xj(e))*(4*x-3*xj(e)-xj(e+1))+...
        c(x)*((2*x-xj(e)-xj(e+1))*(x-xj(e+1)))*((2*x-xj(e+1)-xj(e))*(x-xj(e))),xj(e),xj(e+1));
    K(Index(e,1),Index(e,3))=K(Index(e,1),Index(e,3))+Gauss3(@(x) d(x)*(4*x-3*xj(e+1)-xj(e))*(4*(xj(e+1)-2*x+xj(e)))+c(x)*((2*x-xj(e)-xj(e+1))*(x-xj(e+1)))*(4*(x-xj(e))*(xj(e+1)-x)),xj(e),xj(e+1));
    K(Index(e,3),Index(e,1))=K(Index(e,3),Index(e,1))+Gauss3(@(x) d(x)*(4*x-3*xj(e+1)-xj(e))*(4*(xj(e+1)-2*x+xj(e)))+c(x)*((2*x-xj(e)-xj(e+1))*(x-xj(e+1)))*(4*(x-xj(e))*(xj(e+1)-x)),xj(e),xj(e+1));
    K(Index(e,2),Index(e,3))=K(Index(e,2),Index(e,3))+Gauss3(@(x) d(x)*(4*x-3*xj(e)-xj(e+1))*(4*(xj(e+1)-2*x+xj(e)))+c(x)*((2*x-xj(e+1)-xj(e))*(x-xj(e)))*(4*(x-xj(e))*(xj(e+1)-x)),xj(e),xj(e+1));
    K(Index(e,3),Index(e,2))=K(Index(e,3),Index(e,2))+Gauss3(@(x) d(x)*(4*x-3*xj(e)-xj(e+1))*(4*(xj(e+1)-2*x+xj(e)))+c(x)*((2*x-xj(e+1)-xj(e))*(x-xj(e)))*(4*(x-xj(e))*(xj(e+1)-x)),xj(e),xj(e+1));
end
K(1,1)=K(1,1)+Gauss3(@(x) d(x)*(4*x-3*xj(1)-xj(2))^2+c(x)*((2*x-xj(2)-xj(1))*(x-xj(1)))^2,xj(1),xj(2));
K(1,N)=K(1,N)+Gauss3(@(x) d(x)*(4*x-3*xj(1)-xj(2))*(4*(xj(2)-2*x+xj(1)))+c(x)*((2*x-xj(1)-xj(2))*(x-xj(1)))*(4*(x-xj(1))*(xj(2)-x)),xj(1),xj(2));
K(N,1)=K(N,1)+Gauss3(@(x) d(x)*(4*x-3*xj(1)-xj(2))*(4*(xj(2)-2*x+xj(1)))+c(x)*((2*x-xj(1)-xj(2))*(x-xj(1)))*(4*(x-xj(1))*(xj(2)-x)),xj(1),xj(2));
K(N,N)=K(N,N)+Gauss3(@(x) d(x)*(4*(xj(2)-2*x+xj(1)))^2+c(x)*(4*(x-xj(1))*(xj(2)-x))^2,xj(1),xj(2));
K(N-1,N-1)=K(N-1,N-1)+Gauss3(@(x) d(x)*(4*x-3*xj(N+1)-xj(N))^2+c(x)*((2*x-xj(N)-xj(N+1))*(x-xj(N+1)))^2,xj(N),xj(N+1));
K(N-1,2*N-1)=K(N-1,2*N-1)+Gauss3(@(x) d(x)*(4*x-3*xj(N+1)-xj(N))*(4*(xj(N+1)-2*x+xj(N)))+c(x)*((2*x-xj(N+1)-xj(N))*(x-xj(N+1)))*(4*(x-xj(N))*(xj(N+1)-x)),xj(N),xj(N+1));
K(2*N-1,N-1)=K(2*N-1,N-1)+Gauss3(@(x) d(x)*(4*x-3*xj(N+1)-xj(N))*(4*(xj(N+1)-2*x+xj(N)))+c(x)*((2*x-xj(N+1)-xj(N))*(x-xj(N+1)))*(4*(x-xj(N))*(xj(N+1)-x)),xj(N),xj(N+1));
K(2*N-1,2*N-1)=K(2*N-1,2*N-1)+Gauss3(@(x) d(x)*(4*(xj(N+1)-2*x+xj(N)))^2+c(x)*(4*(x-xj(N))*(xj(N+1)-x))^2,xj(N),xj(N+1));
K=K./h^4;
%Calculate F
F=zeros(2*N-1,1);
for e=2:N-1
    F(Index(e,1))=F(Index(e,1))+Gauss3(@(x) f(x)*((2*x-xj(e)-xj(e+1))*(x-xj(e+1)))/h^2,xj(e),xj(e+1));
    F(Index(e,2))=F(Index(e,2))+Gauss3(@(x) f(x)*((2*x-xj(e+1)-xj(e))*(x-xj(e)))/h^2,xj(e),xj(e+1));
    F(Index(e,3))=F(Index(e,3))+Gauss3(@(x) f(x)*(4*(x-xj(e))*(xj(e+1)-x))/h^2,xj(e),xj(e+1));
end
F(1)=F(1)+Gauss3(@(x) f(x)*((2*x-xj(2)-xj(1))*(x-xj(1)))/h^2,xj(1),xj(2));
F(N)=F(N)+Gauss3(@(x) f(x)*(4*(x-xj(1))*(xj(2)-x))/h^2,xj(1),xj(2));
F(N-1)=F(N-1)+Gauss3(@(x) f(x)*((2*x-xj(N)-xj(N+1))*(x-xj(N+1)))/h^2,xj(N),xj(N+1));
F(2*N-1)=F(2*N-1)+Gauss3(@(x) f(x)*(4*(x-xj(N))*(xj(N+1)-x))/h^2,xj(N),xj(N+1));
%solve the equation
U=K\F;
end