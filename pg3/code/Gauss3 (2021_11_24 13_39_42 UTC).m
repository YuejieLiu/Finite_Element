function value=Gauss3(f,a,b)
value=0;
x=[-sqrt(3/5),0,sqrt(3/5)];
A=[5/9,8/9,5/9];
for k=1:3
    value=value+A(k)*f(((b-a)*x(k)+a+b)/2);
end
value=(b-a)/2*value;
end