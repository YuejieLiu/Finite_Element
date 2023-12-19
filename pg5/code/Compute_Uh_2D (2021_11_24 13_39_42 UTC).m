function Uh=Compute_Uh_2D(x,y,U,element,node,Coeff,C)
[N,~]=size(element);
for m=1:N
    P=element(m,2:4);
    V=node(P+1,2:3);
    [in,on] = inpolygon(x,y,V(:,1),V(:,2));
    if in||on
        break;
    end
end
vertex=element(m,2:4);
a=[];
if ~isempty(find(C==vertex(1),1))
    a=[a,find(C==vertex(1))];
end
if ~isempty(find(C==vertex(2),1))
    a=[a,find(C==vertex(2))];
end
if ~isempty(find(C==vertex(3),1))
    a=[a,find(C==vertex(3))];
end
for k=1:3
    base=[];
    for i=1:numel(a)
        j=find(vertex==C(a(i)));
        base=[base,Coeff(3*m+j-3,1)*x+Coeff(3*m+j-3,2)*y+Coeff(3*m+j-3,3)];
    end
    Uh=dot(U(a),base);
end
end