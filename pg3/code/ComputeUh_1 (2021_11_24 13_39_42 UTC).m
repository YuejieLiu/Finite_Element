function y=ComputeUh_1(U,x)
N=numel(U)+1;
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
    y=U(1)*(x-0)/h;
    return
else
    if j==N-1
        y=U(N-1)*(1-x)/h;
        return
    else
        y=U(j)*((j+1)*h-x)/h+U(j+1)*(x-j*h)/h;
        return
    end
end
end