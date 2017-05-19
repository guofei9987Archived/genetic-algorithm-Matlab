function fitnessY=my_tsp_fitness(Chrom,citycoor,citydist)

%x           优化对象
%citycoor    坐标
%cityDist    距离
[m,n]=size(Chrom);
fitnessY=zeros(m,1);
for i=1:m
    temp=0;
    for j=1:n-1
        temp=temp+citydist(Chrom(i,j),Chrom(i,j+1));
    end
    temp=temp+citydist(Chrom(i,n),Chrom(i,1));
    fitnessY(i,1)=temp;
end
end
