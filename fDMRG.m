function [A,E] = fDMRG(MP,MO,N)
%MP:Initial MPS
% MO: matrix product operator
% N: number of sweeps
%< Output >
% A1 : [cell] MPS of the ground state
% E : [numeric] Ground-state energy of the  Hamiltonian.
%


t = size(MP,2);
HeisMPO = MO;
MP = rightcanon(MP);
 
    
    while(N~=0)                           
        for i= (1:t-1)                        
            if (i==1)
                
                
                Beff = calcR(MP(2:t),HeisMPO(2:t));
                [A1,E]=applyMPOsite([],HeisMPO{1},Beff);
                A1 = reshape(A1, [size(A1,1)*size(A1,2), size(A1,3)]); % Convert the rank 3 tensor into a matrix to perform SVD
                [U,S,V] = svd(A1,'econ'); % econ helps in memory 
                MP{1} = reshape(U,[size(U,1)/size(MP{i},2), size(MP{i},2), size(U,2)]);
                MP{i+1} = ncon({(S*V'),MP{i+1}},{[-1 1],[1 -2 -3]});
                
                
                

                
                
                
            else
                
                Aeff = calcL(MP(1:i-1),HeisMPO(1:i-1));
                Beff = calcR(MP(i+1:t),HeisMPO(i+1:t));
                [A1,E]=applyMPOsite(Aeff,HeisMPO{i},Beff);
                A1 = reshape(A1, [size(A1,1)*size(A1,2), size(A1,3)]); % Convert the rank 3 tensor into a matrix to perform SVD
                [U,S,V] = svd(A1,'econ'); % econ helps in memory 
                MP{i} = reshape(U,[size(U,1)/size(MP{i},2), size(MP{i},2), size(U,2)]);
                MP{i+1} = ncon({(S*V'),MP{i+1}},{[-1 1],[1 -2 -3]});
                
                                                    
            end             
            
        end
        
        % end for loop
          
 MP = leftcanon(MP);       
          
          for i= (t:-1:2)                        
            if (i==t)
                
                
                Aeff = calcL(MP(1:t-1),HeisMPO(1:t-1));
                [A1,E]=applyMPOsite(Aeff,HeisMPO{i},[]);
                
                A1 = reshape(A1, [size(A1,1), size(A1,2)*size(A1,3)]);
                [U,S,V] = svd(A1,'econ');
                MP{i} = reshape(V',[size(V,2), size(MP{i},2), size(V,1)/size(MP{i},2)]);
                MP{i-1} = ncon({(MP{i-1}),U*S},{[-1 -2 1],[1 -3]});
                
                

                
                
                
            else
                
                Aeff = calcL(MP(1:i-1),HeisMPO(1:i-1));
                Beff = calcR(MP(i+1:t),HeisMPO(i+1:t));
                [A1,E]=applyMPOsite(Aeff,HeisMPO{i},Beff);
                A1 = reshape(A1, [size(A1,1), size(A1,2)*size(A1,3)]);
                [U,S,V] = svd(A1,'econ');
                MP{i} = reshape(V',[size(V,2), size(MP{i},2), size(V,1)/size(MP{i},2)]);
                MP{i-1} = ncon({(MP{i-1}),U*S},{[-1 -2 1],[1 -3]});
                
                                                    
            end             
            
         end
          N=N-1;
          %x =E1;
         %[MP,~]=fDMRG(MP,MO,N);
         %[MP,~]=fDMRG(MP,MO,N,x);
    
    
        
       

    
    
    
    end

 A=MP;
 
 %E=x;

end