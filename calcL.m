%Function for evaluating <ps|H|ps>, helper function for the main dmrg
%Ps = cell(1:N) where each element is a tensor 
%outputs a rank 3 tensor 
function [A] = calcL(MP,MO)
if (size(MP,2)~= size(MO,2))
    error('ERR: MPS and MPO are not compatible');
end
N = size(MP,2);
A1= ncon({MP{1},MO{1},conj(MP{1})},{[-1,1,-2],[-3,1,-4,2],[-5,2,-6]});
A1 = reshape(A1,size(A1,6),size(A1,4),size(A1,2));
if (N==1)
        
     A=A1;
else
for i=(2:N)
    X = conj(MP{i});
    
    
   
    
    A1 =ncon({A1,MP{i}},{[1,-1,-2],[1,-3,-4]});
    A1 = ncon({A1,X},{[-1,1,-2,-3],[1,-4,-5]});
    A1 = ncon({MO{i},A1},{[1,2,-2,3],[1,2,-1,3,-3]});
        
        
    
end
end
    
%A = Aeff;
A=A1;
end
