function [B] = calcR(MP,MO)
if (size(MP,2)~= size(MO,2))
    error('ERR: MPS and MPO are not compatible');
end
N = size(MP,2);
Beff= ncon({MP{N},MO{N},conj(MP{N})},{[-1,1,-2],[-3,1,-4,2],[-5,2,-6]});
Beff = reshape(Beff,size(Beff,1),size(Beff,3),size(Beff,5));
if (N==1)
        
     B=Beff;
else
for i=(N-1:-1:1)
    X = conj(MP{i});
    
     Beff= ncon({MP{i},X,MO{i},Beff},{[-1,1,2],[-3,3,4],[-2,1,5,3],[2,5,4]},[2 4 5 1 3]);
    % Beff= ncon({MO{i},Beff},{[-1,-2,1,-3],[-4,1,-5]}); 

     %Beff= ncon({MP{i},Beff},{[-1,1,2],[-2,-3,1,2,-4]});
     %Beff= ncon({X,Beff},{[-1,1,2],[-3,-2,1,2]});
      
        
    
end
end

%B = Beff;
B=Beff;
end
