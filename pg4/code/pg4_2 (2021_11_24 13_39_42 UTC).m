clear all;
N=10;
eps=10^(-7);
error1=0;
order1=0;
error_inf=0;
order_inf=0;
% u=@(x) -(-exp(1/eps)*x.^2+exp(x./eps)+2*eps*(-exp(1/eps)*x+exp(x./eps)+x-1)+x.^2-1)./(2*(exp(1/eps)-1));
b=-1/eps;
u=@(x) x.^2/2+eps*x-(1+2*eps)/2*(exp(x/eps+b)-exp(b))./(exp(b+1/eps)-exp(b));
f=@(x) x;
while N<=320
    tau=1-2*eps*log(N);
    %compute coefficient
    h1=tau/N;
    h2=(1-tau)/N;
    xj=zeros(2*N+1,1);
    for j=1:N+1
        xj(j)=(j-1)*h1;
    end
    for j=N+2:2*N+1
        xj(j)=tau+(j-N-1)*h2;
    end
    Index=zeros(2*N,2);
    for e=2:2*N-1
        Index(e,1)=e-1;
        Index(e,2)=e;
    end
    Index(1,1)=0;
    Index(1,2)=1;
    Index(2*N,1)=2*N-1;
    Index(2*N,2)=0;
    %Assemble
    K=zeros(2*N-1,2*N-1);
    for e=2:N
        K(Index(e,1),Index(e,1))=K(Index(e,1),Index(e,1))+Gauss3(@(x) eps-(xj(e+1)-x),xj(e),xj(e+1));
        K(Index(e,2),Index(e,2))=K(Index(e,2),Index(e,2))+Gauss3(@(x) eps+(x-xj(e)),xj(e),xj(e+1));
        K(Index(e,2),Index(e,1))=K(Index(e,2),Index(e,1))+Gauss3(@(x) -eps-(x-xj(e)),xj(e),xj(e+1));
        K(Index(e,1),Index(e,2))=K(Index(e,1),Index(e,2))+Gauss3(@(x) -eps+(xj(e+1)-x),xj(e),xj(e+1));
    end
    K(1,1)=K(1,1)+Gauss3(@(x) eps+(x-xj(1)),xj(1),xj(2));
    K=K./h1^2;
    K1=zeros(2*N-1,2*N-1);
    for e=N+1:2*N-1
        K1(Index(e,1),Index(e,1))=K1(Index(e,1),Index(e,1))+Gauss3(@(x) eps-(xj(e+1)-x),xj(e),xj(e+1));
        K1(Index(e,2),Index(e,2))=K1(Index(e,2),Index(e,2))+Gauss3(@(x) eps+(x-xj(e)),xj(e),xj(e+1));
        K1(Index(e,2),Index(e,1))=K1(Index(e,2),Index(e,1))+Gauss3(@(x) -eps-(x-xj(e)),xj(e),xj(e+1));
        K1(Index(e,1),Index(e,2))=K1(Index(e,1),Index(e,2))+Gauss3(@(x) -eps+(xj(e+1)-x),xj(e),xj(e+1));
    end
    K1(2*N-1,2*N-1)=K1(2*N-1,2*N-1)+Gauss3(@(x) eps-(xj(2*N+1)-x),xj(2*N),xj(2*N+1));
    K1=K1./h2^2;
    K=K+K1;
    %Calculate F
    F=zeros(2*N-1,1);
    for e=2:N
        F(Index(e,1))=F(Index(e,1))+Gauss3(@(x) f(x)*(xj(e+1)-x)/h1,xj(e),xj(e+1));
        F(Index(e,2))=F(Index(e,2))+Gauss3(@(x) f(x)*(x-xj(e))/h1,xj(e),xj(e+1));
    end
    F(1)=F(1)+Gauss3(@(x) f(x)*(x-xj(1))/h1,xj(1),xj(2));
    for e=N+1:2*N-1
        F(Index(e,1))=F(Index(e,1))+Gauss3(@(x) f(x)*(xj(e+1)-x)/h2,xj(e),xj(e+1));
        F(Index(e,2))=F(Index(e,2))+Gauss3(@(x) f(x)*(x-xj(e))/h2,xj(e),xj(e+1));
    end
    F(2*N-1)=F(2*N-1)+Gauss3(@(x) f(x)*(xj(2*N+1)-x)/h2,xj(2*N),xj(2*N+1));
    %solve the equation
    U=K\F;
    %compute errors and orders
    point=[-sqrt(3/5),0,sqrt(3/5)];
    X=zeros(2*N,3);
    for k=1:N
        X(k,:)=(h1.*point+(2*k-1)*h1)./2;
    end
    for k=N+1:2*N
        X(k,:)=(h2.*point+2*tau+(2*(k-N)-1)*h2)./2;
    end
    U0=u(X);
    Uh=zeros(size(X));
    for i=1:2*N
        for j=1:3
            Uh(i,j)=ComputeUh_Shinshkin(U,X(i,j),tau);
        end
    end
    A=[5/9,8/9,5/9];
    value=0;
    E=abs(U0-Uh);
    for k=1:N
    value=value+dot(A,E(k,:));
    end
    value=h1/2*value;
    val=0;
    for k=N+1:2*N
        val=val+dot(A,E(k,:));
    end
    val=h2/2*val;
    value=value+val;
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