N=10;
error1=0;
order1=0;
error_inf=0;
order_inf=0;
while N<=80
    h=1/N;
    %K is the stiffness matrix
    %We use Simpson rule to compute all of the following integrations
    K=zeros(2*N-1,2*N-1);
    for j=1:N-2
        K(j,j)=14/(3*h);
        K(j,j+1)=1/(3*h);
        K(j+1,j)=1/(3*h);
        K(j,j+N-1)=-8/(3*h);
        K(j,j+N)=-8/(3*h);
        K(j+N-1,j)=-8/(3*h);
        K(j+N,j)=-8/(3*h);
    end
    K(N-1,N-1)=14/(3*h);
    K(N-1,N-1+N-1)=-8/(3*h);
    K(N-1,N-1+N)=-8/(3*h);
    K(N-1+N-1,N-1)=-8/(3*h);
    K(N-1+N,N-1)=-8/(3*h);
    for j=N:2*N-1
        K(j,j)=16/(3*h);
    end
    %U is the vector of function value on grid points
    %F is the inner product of f and basis function phi_i
    %-u''=f=-2cosx+(x-1)sinx
    f=@(x) -2*cos(x)+(x-1)*sin(x);
    F=zeros(2*N-1,1);
    for j=1:N-1
        F(j)=h*f(j*h)/3;
    end
    for j=N:2*N-1
        F(j)=2/3*h*f((2*j-2*N+1)*h/2);
    end
    %solve linear equations to get the coefficcient of u_h
    U=K\F;
    %compute errors and orders
    X=1/320*(1:320);
    U0=(X-1).*sin(X);
    Uh=zeros(size(X));
    for j=1:320
        Uh(j)=ComputeUh_2(U,X(j));
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
    fprintf('N=%d,L1 error=%4.12f,order1=%4.12f,L_inf error=%4.12f,order_inf=%4.12f \n',N, error_temp_1,order1,...
        error_temp_inf,order_inf);
    N=N*2;
end

