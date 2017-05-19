%二进制转实值
%!需要改进
%飞飞出品QQ513829987
%Chrom 二进制矩阵，表示一个种群
%FieldD=[len;lb;ub;code;scale;lbin;ubin]
%sum(len)=length(Chrom)每个子串的长度
%lb每个变量的下界，ub每个变量的上界
%code(i)=1  二进制编码，code(i)=0  格雷编码
%scale二进制向量，scale(i)=1表示算术刻度，scale(i)=0表示对数刻度
%lbin这个值为0表示范围中去掉下界，为1表示包含下界
%ubin为0表示去掉上界，为1表示包含上界

%
% len = [50,50];%每个子串（一个子串对应一个变量）的长度
% lb =[1,2];%每个变量的下界
% ub = [5,6];%每个变量的上界
% code = [0,0];%code(i)=0表示用格雷编码，code(i)=1表示用二进制
% scale = [0,0];%scale(i)=0算术刻度,scale(i)=1对数刻度
% lin = [1,1];%lin(i)=0表示不含边界，=1表示含边界
% uin = [1,1];%同上
% FieldD=[len
%     lb
%     ub
%     code
%     scale
%     lin
%     uin]
%输出Phen是实数

function Phen = bs2rv(Chrom,FieldD)

[Nind,~] = size(Chrom);


[~,Nvar] = size(FieldD);


len = FieldD(1,:);%每个子串（一个子串对应一个变量）的长度
lb = FieldD(2,:);%每个变量的下界
ub = FieldD(3,:);%每个变量的上界
code = ~(~FieldD(4,:));%code(i)=0表示用格雷编码，code(i)=1表示用二进制
scale = ~(~FieldD(5,:));%scale(i)=0算术刻度,scale(i)=1对数刻度；注意对数刻度基因不能为0
lin = ~(~FieldD(6,:));%lin(i)=0表示不含边界，=1表示含边界
uin = ~(~FieldD(7,:));%同上

Phen = zeros(Nind,Nvar);

lf = cumsum(len);
li = cumsum([1 len]);
Prec = .5 .^ len;

logsgn = sign(lb(scale));

%对数刻度2算术刻度
lb(scale) = log( abs(lb(scale)) );
ub(scale) = log( abs(ub(scale)) );
delta = ub - lb;

num = (~lin) .* Prec;
den = (lin + uin - 1) .* Prec;

for i = 1:Nvar,
    idx = li(i):lf(i);
    if code(i) % Gray decoding
        Chrom(:,idx)=rem(cumsum(Chrom(:,idx)')',2);
    end
    Phen(:,i) = Chrom(:,idx) * [ (.5).^(1:len(i))' ];
    Phen(:,i) = lb(i) + delta(i) * (Phen(:,i) + num(i)) ./ (1 - den(i));
end

expand = ones(Nind,1);
if any(scale)
    Phen(:,scale) = logsgn(expand,:) .* exp(Phen(:,scale));
end
end
