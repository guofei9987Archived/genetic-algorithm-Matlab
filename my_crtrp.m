%生成种群--实数形式
%飞飞出品%QQ513829987
%Chrom:种群size=(pop,lind)
%pop：个体数量
%lind：基因数量
%base：2*lind矩阵,分别表示上界和下界，每个基因的可能取值(2,2,2,2,3,4
%                                                     1,1,1,1,1,0.5)
function Chrom=my_crtrp(pop,lind,base)
if nargin==2
    base_lb=2*ones(1,lind);
    base_ub=zeros(1,lind);
else base_lb=base(1,:);
    base_ub=base(2,:);
end

if size(base_lb,2)~=lind
    error('base的列数就是lind')
end

Chrom=rand(pop,lind).*(ones(pop,1)*base_lb)+ones(pop,1)*base_ub;


end