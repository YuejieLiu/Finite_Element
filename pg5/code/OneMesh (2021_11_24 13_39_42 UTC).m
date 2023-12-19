function [L1,Linf,N]=OneMesh(StrNode,StrEle) 
    num_n=4;
    num_e=13;
    f=@(x,y) 8*pi*pi*sin(2*pi*x).*sin(2*pi*y);
    u=@(x,y) sin(2*pi*x).*sin(2*pi*y);
    %importdata是按行读取，结果是向量
    node=importdata(StrNode);
    element=importdata(StrEle);
    N=element(1);
    %index of node starts from 0
    node=reshape(node(2:end),num_n,[])';
    element=reshape(element(2:end),num_e,[])';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %compute basis function
    Coeff=zeros(3*N,3);
    Area=zeros(N,1);
    for m=1:N
        NodeIndex=element(m,2:4);
        b1=[1,0,0]';
        b2=[0,1,0]';
        b3=[0,0,1]';
        A=node(NodeIndex+1,2:3);
        A(:,3)=[1,1,1]';
        Coeff(3*m-2,:)=A\b1;
        Coeff(3*m-1,:)=A\b2;
        Coeff(3*m,:)=A\b3;
        Area(m)=abs(det(A))/2;
    end
    %compute the number of inner points
    B=FindBoundary(element,N);
    C=setdiff(node(:,1),B);%index of inner points/base
    M=numel(C);
    %Assemble stiffness matrix
    K=zeros(M,M);
    for m=1:N
        %若一条边在边界上，只计算对点nodal basis 的积分
        if element(m,5)==-1
            i=find(C==element(m,2));
            K(i,i)=K(i,i)+(Coeff(3*m-2,1)^2+Coeff(3*m-2,2)^2)*Area(m);
        else
            if element(m,6)==-1
                i=find(C==element(m,3));
                K(i,i)=K(i,i)+(Coeff(3*m-1,1)^2+Coeff(3*m-1,2)^2)*Area(m);
            else
                if element(m,7)==-1
                    i=find(C==element(m,4));
                    K(i,i)=K(i,i)+(Coeff(3*m,1)^2+Coeff(3*m,2)^2)*Area(m);
                else%没有边在边界上
                    %如果有一个点在边界上
                    if ismember(element(m,2),B)
                        i=find(C==element(m,3));
                        j=find(C==element(m,4));
                        K(i,i)=K(i,i)+(Coeff(3*m-1,1)^2+Coeff(3*m-1,2)^2)*Area(m);
                        K(j,j)=K(j,j)+(Coeff(3*m,1)^2+Coeff(3*m,2)^2)*Area(m);
                        K(i,j)=K(i,j)+(Coeff(3*m-1,1)*Coeff(3*m,1)+Coeff(3*m-1,2)*Coeff(3*m,2))*Area(m);
                        K(j,i)=K(j,i)+(Coeff(3*m-1,1)*Coeff(3*m,1)+Coeff(3*m-1,2)*Coeff(3*m,2))*Area(m);
                    else
                        if ismember(element(m,3),B)
                            i=find(C==element(m,2));
                            j=find(C==element(m,4));
                            K(i,i)=K(i,i)+(Coeff(3*m-2,1)^2+Coeff(3*m-2,2)^2)*Area(m);
                            K(j,j)=K(j,j)+(Coeff(3*m,1)^2+Coeff(3*m,2)^2)*Area(m);
                            K(i,j)=K(i,j)+(Coeff(3*m-2,1)*Coeff(3*m,1)+Coeff(3*m-2,2)*Coeff(3*m,2))*Area(m);
                            K(j,i)=K(j,i)+(Coeff(3*m-2,1)*Coeff(3*m,1)+Coeff(3*m-2,2)*Coeff(3*m,2))*Area(m);
                        else
                            if ismember(element(m,4),B)
                                i=find(C==element(m,2));
                                j=find(C==element(m,3));
                                K(i,i)=K(i,i)+(Coeff(3*m-2,1)^2+Coeff(3*m-2,2)^2)*Area(m);
                                K(j,j)=K(j,j)+(Coeff(3*m-1,1)^2+Coeff(3*m-1,2)^2)*Area(m);
                                K(i,j)=K(i,j)+(Coeff(3*m-2,1)*Coeff(3*m-1,1)+Coeff(3*m-2,2)*Coeff(3*m-1,2))*Area(m);
                                K(j,i)=K(j,i)+(Coeff(3*m-2,1)*Coeff(3*m-1,1)+Coeff(3*m-2,2)*Coeff(3*m-1,2))*Area(m);
                            else
                                i=find(C==element(m,2));
                                j=find(C==element(m,3));
                                k=find(C==element(m,4));
                                K(i,i)=K(i,i)+(Coeff(3*m-2,1)^2+Coeff(3*m-2,2)^2)*Area(m);
                                K(j,j)=K(j,j)+(Coeff(3*m-1,1)^2+Coeff(3*m-1,2)^2)*Area(m);
                                K(k,k)=K(k,k)+(Coeff(3*m,1)^2+Coeff(3*m,2)^2)*Area(m);
                                K(i,j)=K(i,j)+(Coeff(3*m-2,1)*Coeff(3*m-1,1)+Coeff(3*m-2,2)*Coeff(3*m-1,2))*Area(m);
                                K(j,i)=K(j,i)+(Coeff(3*m-2,1)*Coeff(3*m-1,1)+Coeff(3*m-2,2)*Coeff(3*m-1,2))*Area(m);
                                K(i,k)=K(i,k)+(Coeff(3*m-2,1)*Coeff(3*m,1)+Coeff(3*m-2,2)*Coeff(3*m,2))*Area(m);
                                K(k,i)=K(k,i)+(Coeff(3*m-2,1)*Coeff(3*m,1)+Coeff(3*m-2,2)*Coeff(3*m,2))*Area(m);
                                K(j,k)=K(j,k)+(Coeff(3*m-1,1)*Coeff(3*m,1)+Coeff(3*m-1,2)*Coeff(3*m,2))*Area(m);
                                K(k,j)=K(k,j)+(Coeff(3*m-1,1)*Coeff(3*m,1)+Coeff(3*m-1,2)*Coeff(3*m,2))*Area(m);
                            end
                        end
                    end     
                end
            end
        end
    end
    %compute F
    F=zeros(M,1);
    for m=1:N
        NodeIndex=element(m,2:4);
        A=node(NodeIndex+1,2:3);
        v1=@(x,y) Coeff(3*m-2,1)*x+Coeff(3*m-2,2)*y+Coeff(3*m-2,3);
        v2=@(x,y) Coeff(3*m-1,1)*x+Coeff(3*m-1,2)*y+Coeff(3*m-1,3);
        v3=@(x,y) Coeff(3*m,1)*x+Coeff(3*m,2)*y+Coeff(3*m,3);
        %若一条边在边界上，只计算对点nodal basis 的积分
        if element(m,5)==-1
            i=find(C==element(m,2));
            F(i)=F(i)+Gauss2D(@(x,y) f(x,y)*v1(x,y),A(1,1),A(1,2),A(2,1),A(2,2),A(3,1),A(3,2));
        else
            if element(m,6)==-1
                i=find(C==element(m,3));
                F(i)=F(i)+Gauss2D(@(x,y) f(x,y)*v2(x,y),A(1,1),A(1,2),A(2,1),A(2,2),A(3,1),A(3,2));
            else
                if element(m,7)==-1
                    i=find(C==element(m,4));
                    F(i)=F(i)+Gauss2D(@(x,y) f(x,y)*v3(x,y),A(1,1),A(1,2),A(2,1),A(2,2),A(3,1),A(3,2));
                else%没有边在边界上
                    %如果有一个点在边界上
                    if ismember(element(m,2),B)
                        i=find(C==element(m,3));
                        j=find(C==element(m,4));
                        F(i)=F(i)+Gauss2D(@(x,y) f(x,y)*v2(x,y),A(1,1),A(1,2),A(2,1),A(2,2),A(3,1),A(3,2));
                        F(j)=F(j)+Gauss2D(@(x,y) f(x,y)*v3(x,y),A(1,1),A(1,2),A(2,1),A(2,2),A(3,1),A(3,2));
                    else
                        if ismember(element(m,3),B)
                            i=find(C==element(m,2));
                            j=find(C==element(m,4));
                            F(i)=F(i)+Gauss2D(@(x,y) f(x,y)*v1(x,y),A(1,1),A(1,2),A(2,1),A(2,2),A(3,1),A(3,2));
                            F(j)=F(j)+Gauss2D(@(x,y) f(x,y)*v3(x,y),A(1,1),A(1,2),A(2,1),A(2,2),A(3,1),A(3,2));
                        else
                            if ismember(element(m,4),B)
                                i=find(C==element(m,2));
                                j=find(C==element(m,3));
                                F(i)=F(i)+Gauss2D(@(x,y) f(x,y)*v1(x,y),A(1,1),A(1,2),A(2,1),A(2,2),A(3,1),A(3,2));
                                F(j)=F(j)+Gauss2D(@(x,y) f(x,y)*v2(x,y),A(1,1),A(1,2),A(2,1),A(2,2),A(3,1),A(3,2));
                            else
                                i=find(C==element(m,2));
                                j=find(C==element(m,3));
                                k=find(C==element(m,4));
                                F(i)=F(i)+Gauss2D(@(x,y) f(x,y)*v1(x,y),A(1,1),A(1,2),A(2,1),A(2,2),A(3,1),A(3,2));
                                F(j)=F(j)+Gauss2D(@(x,y) f(x,y)*v2(x,y),A(1,1),A(1,2),A(2,1),A(2,2),A(3,1),A(3,2));
                                F(k)=F(k)+Gauss2D(@(x,y) f(x,y)*v3(x,y),A(1,1),A(1,2),A(2,1),A(2,2),A(3,1),A(3,2));
                            end
                        end
                    end
                end
            end
        end
    end
    %solve coefficient
    U=K\F;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Generate test points (Our code here is redundant, because we repeatly
    %compute it in Gauss2D, which could be improved at some time)
    [L1,Linf,img]=ComputeError(U,element,node,Coeff,C,u);
    figure
    plot3(img(:,1),img(:,2),img(:,3),'.');
    title([StrNode,StrEle,'Error']);
    %尽管我们用高斯结点计算误差，但是我们仍然用均匀网格画出u和uh的图像
    [x,y]=meshgrid(0:0.01:1,0:0.01:1);
    [row,col]=size(x);
    z=zeros(row,col);
    for i=1:row
        for j=1:col
            z(i,j)=Compute_Uh_2D(x(i,j),y(i,j),U,element,node,Coeff,C);
        end
    end
    figure
    mesh(x,y,z);
    title([StrNode,StrEle,'FEM Solution']);
end