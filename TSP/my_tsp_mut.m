function A=my_tsp_mut(A,n1,n2)
temp=A(1,n1);
A(1,n1)=A(1,n2);
A(1,n2)=temp;
end