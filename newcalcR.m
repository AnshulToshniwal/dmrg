




function [B] = newcalcR(MP,MO)
if (size(MP,2)~= size(MO,2))
    error('ERR: MPS and MPO are not compatible');
end
N = size(MP,2);
Beff = cell(1,N);
for i=(1:N)
    X = conj(MP{i});
    
    % evaluating the expression for individual sites and replacing with a
    % rank 6 tensor
  %  Beff{i}= ncon({X,MO{i},MP{i}},{[-1,1,-2],[-3,1,-4,2],[-5,2,-6]});
    Beff{i}= ncon({MP{i},MO{i},X},{[-1,1,-2],[-3,1,-4,2],[-5,2,-6]});
    
        
        
    
end
if (N==1)
        
     Z = Beff{1};
     Z = reshape(Z,size(Z,1),size(Z,3),size(Z,5));
else
    for i=(1:N-1)
        if (i==1)
            Z = ncon({Beff{i},Beff{i+1}},{[-1,1,-2,2,-3,3],[1,-4,2,-5,3,-6]});
           
        else
            Z = ncon({Z,Beff{i+1}},{[-1,-2,-3,1,2,3],[1,-4,2,-5,3,-6]});
            
        end    
    end
    Z = reshape(Z,size(Z,1),size(Z,2),size(Z,3));
end

%B = Beff;
B=Z;
end