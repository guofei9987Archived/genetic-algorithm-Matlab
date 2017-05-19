%生成种群-离散形式
%飞飞出品QQ513829987
%Chrom:种群size=(pop,lind)
%pop：个体数量
%lind：基因数量
%base：1*lind矩阵，每个基因的可能取值(2,2,2,2,3,4)
function Chrom=my_crtbp(pop,lind,base)
if nargin==2
    base=2*ones(1,lind);
end

if size(base,2)~=lind
    error('base的列数就是lind')
end

Chrom=floor(rand(pop,lind).*(ones(pop,1)*base));

end