function [A,B]=my_tsp_crossover(A,B,n1,n2)
for i=n1:n2
    x=find(A==B(1,i));
    y=find(B==A(1,i));
    [A(1,i),B(1,i)]=exchange(A(1,i),B(1,i));
    [A(1,x),B(1,y)]=exchange(A(1,x),B(1,y));
end

end

%¶Ôµ÷
function [x,y]=exchange(x,y)
temp=x;
x=y;
y=temp;
end