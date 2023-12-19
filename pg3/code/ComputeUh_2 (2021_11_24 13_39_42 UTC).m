function y=ComputeUh_2(U,x)
N=(numel(U)+1)/2;
h=1/N;
j=0;
%find which part x locates at
for i=0:N-1
    if(x>=i*h&&x<=(i+1)*h)
        j=i;
        break;
    end
end
if j==0
    y=U(1)*(2*x-j*h-(j+1)*h)*(x-j*h)/h^2+U(N)*4*(x-0)*(h-x)/h^2;
    return
else
    if j==N-1
        y=U(N-1)*(2*x-j*h-(j+1)*h)*(x-(j+1)*h)/h^2+U(2*N-1)*4*(x-j*h)*(N*h-x)/h^2;
        return
    else
        y=U(j)*(2*x-j*h-(j+1)*h)*(x-(j+1)*h)/h^2+U(j+1)*(2*x-j*h-(j+1)*h)*(x-j*h)/h^2+U(j+N)*4*(x-j*h)*((j+1)*h-x)/h^2;
        return
    end
end
end