%%
%生成问题
clear
clc
close all
num_city=18
citycoor=random('unif',0,1,num_city,2);%生成城市坐标
citydist=squareform(pdist(citycoor));%计算距离矩阵



%%
%算法实现
pop_size=300;%种群
lind=num_city;%基因长度

[~,Chrom]=sort(random('unif',0,1,pop_size,lind),2);


for i=1:500
    
    fitnessY=my_tsp_fitness(Chrom,citycoor,citydist);

    g=ranking(fitnessY);
    %scaling
    
    NewChrIx=rws(g,pop_size);%新个体的索引值
    Chrom=Chrom(NewChrIx,:);%新个体
    %reins
    for j=1:2:pop_size-1
        if rand<0.4
            cc=random('unid',lind,1,2);
            c1=cc(1);
            c2=cc(2);
            [Chrom(j,:),Chrom(j+1,:)]=my_tsp_crossover(Chrom(j,:),Chrom(j+1,:),c1,c2);
        end
    end
    %Chrom=xovmp(Chrom,0.8,1);%交叉
    
    for j=1:pop_size
        if rand<0.2
            cc=random('unid',lind,1,2);
            c1=cc(1);
            c2=cc(2);
            Chrom(j,:)=my_tsp_mut(Chrom(j,:),c1,c2);
        end
    end
    %变异
    
    
    
    
    
end
fitnessY=my_tsp_fitness(Chrom,citycoor,citydist);
[fitnessY_bset,ind]=min(fitnessY);
x=Chrom(ind,:)

figure 

plot(citycoor(:,1),citycoor(:,2),'o')
hold on
plot(citycoor([x,x(1)],1),citycoor([x,x(1)],2),'-')
axis([0,1,0,1])