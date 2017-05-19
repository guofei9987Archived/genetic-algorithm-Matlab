%函数最优问题
clear;clc;
%%
%计算必须的精度
a=0%区间左
b=10%区间右
k=4%精度

con=0;lind=0;
while con==0
    if 2^lind>(b-a)*10^k
        con=1;
    else lind=lind+1;
    end
end
lind

%lind是染色体长度


%%

pop=40
lind=50;

BaseV=2*ones(1,lind);
Chrom=my_crtbp(pop,lind,BaseV);

len = [25,25];%每个子串（一个子串对应一个变量）的长度
lb =[-2,-2];%每个变量的下界
ub = [2,2];%每个变量的上界
code = [0,0];%code(i)=0表示用格雷编码，code(i)=1表示用二进制
scale = [0,0];%scale(i)=0算术刻度,scale(i)=1对数刻度
lin = [1,1];%lin(i)=0表示不含边界，=1表示含边界
uin = [1,1];%同上
FieldD=[len
    lb
    ub
    code
    scale
    lin
    uin];

for i=1:200
Phen=bs2rv(Chrom,FieldD);

fitnessY=myfun(Phen(:,1),Phen(:,2));

g=ranking(fitnessY);
%scaling

NewChrIx=rws(g,pop);%新个体的索引值
Chrom=Chrom(NewChrIx,:);%新个体
%reins

Chrom=xovmp(Chrom,0.8,1);%交叉

Chrom=mut(Chrom,0.01,BaseV);
end

%结果输出
Phen=bs2rv(Chrom,FieldD);
fitnessY=myfun(Phen(:,1),Phen(:,2));
[fitnessY_bset,ind]=min(fitnessY);
Phen(ind,:)