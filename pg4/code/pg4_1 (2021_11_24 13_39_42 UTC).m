clear all;
N=20;
% eps=10^(-1);
eps=10^(-7);
error1=0;
order1=0;
error_inf=0;
order_inf=0;
% u=@(x) -(-exp(1/eps)*x.^2+exp(x./eps)+2*eps*(-exp(1/eps)*x+exp(x./eps)+x-1)+x.^2-1)./(2*(exp(1/eps)-1));
b=-1/eps;
u=@(x) x.^2/2+eps*x-(1+2*eps)/2*(exp(x/eps+b)-exp(b))./(exp(b+1/eps)-exp(b));
f=@(x) x;
while N<=640
    %compute coefficient
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
        K(Index(e,1),Index(e,1))=K(Index(e,1),Index(e,1))+Gauss3(@(x) eps-(xj(e+1)-x),xj(e),xj(e+1));
        K(Index(e,2),Index(e,2))=K(Index(e,2),Index(e,2))+Gauss3(@(x) eps+(x-xj(e)),xj(e),xj(e+1));
        K(Index(e,2),Index(e,1))=K(Index(e,2),Index(e,1))+Gauss3(@(x) -eps-(x-xj(e)),xj(e),xj(e+1));
        K(Index(e,1),Index(e,2))=K(Index(e,1),Index(e,2))+Gauss3(@(x) -eps+(xj(e+1)-x),xj(e),xj(e+1));
    end
    K(1,1)=K(1,1)+Gauss3(@(x) eps+(x-xj(1)),xj(1),xj(2));
    K(N-1,N-1)=K(N-1,N-1)+Gauss3(@(x) eps-(xj(N+1)-x),xj(N),xj(N+1));
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
    %compute errors and orders
    point=[-sqrt(3/5),0,sqrt(3/5)];
    X=zeros(N,3);
    for k=1:N
        X(k,:)=(h.*point+(2*k-1)*h)./2;
    end
    U0=u(X);
    Uh=zeros(size(X));
    for i=1:N
        for j=1:3
            Uh(i,j)=ComputeUh_1(U,X(i,j));
        end
    end
    A=[5/9,8/9,5/9];
    value=0;
    E=abs(U0-Uh);
    for k=1:N
    value=value+dot(A,E(k,:));
    end
    value=h/2*value;
    %draw comparison graph
    U0=reshape(U0',[],1);
    Uh=reshape(Uh',[],1);
    X=reshape(X',[],1);
    figure
    plot(X,U0,X,Uh,'-.');
    title(['N=',num2str(N)]);
    legend('precise solution','numerical solution','Location','northwest');
    error_temp_1=value;
    error_temp_inf=max(E,[],'all');
    if N~=10
        order1=log2(error1/error_temp_1);
        order_inf=log2(error_inf/error_temp_inf);
    end
    error1=error_temp_1;
    error_inf=error_temp_inf;
    fprintf('N=%d,L1 error=%4.12e,order1=%4.12f,L_inf error=%4.12e,order_inf=%4.12f \n',N, error_temp_1,order1,...
        error_temp_inf,order_inf);
    N=N*2;
end