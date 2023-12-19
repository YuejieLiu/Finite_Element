function y=ComputeUh_Shinshkin(U,x,tau)
N=(numel(U)+1)/2;
h1=tau/N;
h2=(1-tau)/N;
j=-1;
%find which part x locates at
for i=0:N-1
    if(x>=i*h1&&x<=(i+1)*h1)
        j=i;
        break;
    end
end
if j~=-1  %x locates in the first N intervals
    if j==0
        y=U(1)*(x-0)/h1;
        return
    else
        y=U(j)*((j+1)*h1-x)/h1+U(j+1)*(x-j*h1)/h1;
        return
    end
else %x locates in the second N intervals
    for i=N:2*N-1
        if(x>=tau+(i-N)*h2&&x<=tau+(i+1-N)*h2)
            j=i;
            break;
        end
    end
    if j==2*N-1
        y=U(2*N-1)*(1-x)/h2;
        return
    else
        y=U(j)*(tau+(j+1-N)*h2-x)/h2+U(j+1)*(x-(tau+(j-N)*h2))/h2;
        return
    end
end


end