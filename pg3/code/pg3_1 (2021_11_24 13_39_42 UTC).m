N=10;
eps=10^(-1);
%eps=10^(-7);
error1=0;
order1=0;
error_inf=0;
order_inf=0;
u=@(x) x.*(x-1);
f=@(x) -cos(x)*(2*x-1)-2*(sin(x)+2)+c(x).*u(x);
% u=@(x) sin(x).*(x-1);
% f=@(x) -cos(x)*(cos(x)*(x-1)+sin(x))-d(x)*(-sin(x)*(x-1)+2*cos(x))+c(x).*u(x);
while N<=320
    U=ComputeCoef(N,d,c,f);
    %compute errors and orders
    point=[-sqrt(3/5),0,sqrt(3/5)];
    X=zeros(N,3);
    h=1/N;
    for k=1:N
        X(k,:)=(h.*point+(2*k-1)*h)./2;
    end
    U0=X.*(X-1);
    %U0=(X-1).*sin(X);
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
    legend('precise solution','numerical solution');
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