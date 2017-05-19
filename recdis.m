function NewChrom = recdis(OldChrom, XOVR);

% Identify the population size (Nind) and the number of variables (Nvar)
   [Nind,Nvar] = size(OldChrom);

% Identify the number of matings
   Xops = floor(Nind/2);

% which parent gives the value
   Mask1 = (rand(Xops,Nvar)<0.5);
   Mask2 = (rand(Xops,Nvar)<0.5);

% Performs crossover
   odd = 1:2:Nind-1;
   even= 2:2:Nind;
   NewChrom(odd,:)  = (OldChrom(odd,:).* Mask1) + (OldChrom(even,:).*(~Mask1));
   NewChrom(even,:) = (OldChrom(odd,:).* Mask2) + (OldChrom(even,:).*(~Mask2));

% If the number of individuals is odd, the last individual cannot be mated
% but must be included in the new population
   if rem(Nind,2),  NewChrom(Nind,:)=OldChrom(Nind,:); end


% End of function

