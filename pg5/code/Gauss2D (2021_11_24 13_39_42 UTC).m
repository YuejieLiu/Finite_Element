function I=Gauss2D(f,x1,y1,x2,y2,x3,y3)
x=[-sqrt(3/5),0,sqrt(3/5)];
A=[5/9,8/9,5/9];
B=[0,0,1;
    0,1,1;
    1,0,1];
X=[x1,x2,x3]';
Y=[y1,y2,y3]';
X1=B\X;
X2=B\Y;
a1=X1(1);
a2=X1(2);
b1=X1(3);
a3=X2(1);
a4=X2(2);
b2=X2(3);
I=0;
for m=1:3
    for k=1:3
        I=I+A(m)*A(k)*f(a1*(1+x(k))/2+a2*(1-x(k))*(1+x(m))/4+b1,a3*(1+x(k))/2+a4*(1-x(k))*(1+x(m))/4+b2)*(1-x(k))/8;
    end
end
I=I*abs(a1*a4-a2*a3);
return
end