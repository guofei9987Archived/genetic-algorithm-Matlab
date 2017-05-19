function NewChrom = mut(OldChrom,Pm,BaseV)


[Nind, Lind] = size(OldChrom) ;

if nargin < 2, Pm = 0.7/Lind ; end

if (nargin < 3), BaseV = 2*ones(1,Lind);  end

BaseM = BaseV(ones(Nind,1),:);


NewChrom = rem(OldChrom+(rand(Nind,Lind)<Pm).*ceil(rand(Nind,Lind).*(BaseM-1)),BaseM);
