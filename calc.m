function [A] = calc(MP)
N= size(MP,2);
Aeff = cell(1,N);
for i=(1:N)
    X = conj(MP{i});
    Aeff{i}= ncon({MP{i},X},{[-1,1,-2],[-3,1,-4]});
  
        
        
    
end
if (N==1)
        
     Z = Aeff{1};
     
else
    for i=(1:N-1)
        
        if (i==1)
            
            Z = ncon({Aeff{i},Aeff{i+1}},{[-1,1,-2,2],[1,-3,2,-4]});
           % Z = reshape(Z,size(Z,6),size(Z,5),size(Z,4));
        else
            
            
            Z = ncon({Z,Aeff{i+1}},{[-1,-2,1,2],[1,-3,2,-4]}); % some weird error
            %Z = reshape(Z,size(Z,6),size(Z,5),size(Z,4)); % reshape into a rank 3 tensor. The first 3 degrees are 1
    
        
       end
    end
     
end

A=Z;
end