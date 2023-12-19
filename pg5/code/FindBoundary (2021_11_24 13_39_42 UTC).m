function B=FindBoundary(element,N)
%element has the same format as EasyMesh without the first line
%N is the number of element
%B is a vector containing the index of boundary node
B=[];
for m=1:N
    if element(m,5)==-1
        B=[B,element(m,3),element(m,4)];
    else
        if element(m,6)==-1
            B=[B,element(m,2),element(m,4)];
        else
            if element(m,7)==-1
                B=[B,element(m,2),element(m,3)];
            end
        end
    end
end
B=unique(B);
return
end