N=10;
error1=0;
order1=0;
error_inf=0;
order_inf=0;
while N<=80
    h=1/N;
    %K is the stiffness matrix
    K=zeros(N-1,N-1);
    for j=1:N-2
        K(j,j)=2/h;
        K(j,j+1)=-1/h;
        K(j+1,j)=-1/h;
    end
    K(N-1,N-1)=2/h;
    %U is the vector of function value on grid points
    %F is the inner product of f and basis function phi_i
    %-u''=f=-2cosx+(x-1)sinx
    %we use trapozoid rule to calculate integrations
    f=@(x) -2*cos(x)+(x-1)*sin(x);
    F=zeros(N-1,1);
    for j=1:N-1
        F(j)=h*f(j*h);
    end
    %solve linear equations to get the coefficcient of u_h
    U=K\F;
    %compute errors and orders
    X=1/320*(1:320);
    U0=(X-1).*sin(X);
    Uh=zeros(size(X));
    for j=1:320
        Uh(j)=ComputeUh_1(U,X(j));
    end
    %draw comparison graph
    figure
    plot(X,U0,X,Uh,'-.');
    title(['N=',num2str(N)]);
    legend('precise solution','numerical solution');
    error_temp_1=norm(Uh-U0,1);
    error_temp_inf=norm(Uh-U0,Inf);
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

