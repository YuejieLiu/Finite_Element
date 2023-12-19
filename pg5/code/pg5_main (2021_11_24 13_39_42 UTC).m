clear all;
StrNode=char('MESH/gd0.n','MESH/gd1.n','MESH/gd2.n','MESH/gd3.n','MESH/gd4.n');
StrEle=char('MESH/gd0.e','MESH/gd1.e','MESH/gd2.e','MESH/gd3.e','MESH/gd4.e');
[x,y]=meshgrid(0:0.01:1,0:0.01:1);
u=@(x,y) sin(2*pi*x).*sin(2*pi*y);
figure
mesh(x,y,u(x,y));
title('Precise Solution');
error1=0;
order1=0;
error_inf=0;
order_inf=0;
Ek=0;
for p=1:5
    [error_temp_1,error_temp_inf,E_temp]=OneMesh(StrNode(p,:),StrEle(p,:));
    if p~=1
        order1=2*log(error1/error_temp_1)/log(E_temp/Ek);
        order_inf=2*log(error_inf/error_temp_inf)/log(E_temp/Ek);
    end
    error1=error_temp_1;
    error_inf=error_temp_inf;
    Ek=E_temp;
    fprintf('gd%d,L1 error=%4.12e,order1=%4.12f,L_inf error=%4.12e,order_inf=%4.12f \n',p-1, error_temp_1,order1,...
        error_temp_inf,order_inf);
end