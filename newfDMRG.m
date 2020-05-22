function [A,E] = newfDMRG(MP,MO,N,Ls,Rs)
%MP:Initial MPS
% MO: matrix product operator
% N: number of sweeps
%< Output >
% A1 : [cell] MPS of the ground state
% E : [numeric] Ground-state energy of the  Hamiltonian.
% Ls: cell of L tesnors
% Rs: cell of R tensors

Si = size(MP,2);
t = size(MP,2)/2;
HeisMPO = MO;

 
    
    %if(N~=0 && x ~=0)     % The half sweep of the iDMRG result to reach an edge                      
        for i= (t+1:2*t-1)                        
                [A1,E]=applyMPOsite(Ls{i-1},HeisMPO{i},Rs{Si-i});
                A1 = reshape(A1, [size(A1,1)*size(A1,2), size(A1,3)]); % Convert the rank 3 tensor into a matrix to perform SVD
                [U,S,V] = svd(A1,'econ'); % econ helps in memory 
                MP{i} = reshape(U,[size(U,1)/size(MP{i},2), size(MP{i},2), size(U,2)]);
                MP{i+1} = ncon({(S*V'),MP{i+1}},{[-1 1],[1 -2 -3]});
                
                
                if (i<(2*t-1))
                L =ncon({Ls{i-1},MP{i}},{[1,-1,-2],[1,-3,-4]});
                L = ncon({L,conj(MP{i})},{[-1,1,-2,-3],[1,-4,-5]});
                L = ncon({HeisMPO{i},L},{[1,2,-2,3],[1,2,-1,3,-3]});    
                Ls{i} = L;
                R= ncon({MP{i+1},conj(MP{i+1}),HeisMPO{i+1},Rs{Si-i-1}},{[-1,1,2],[-3,3,4],[-2,1,5,3],[2,5,4]},[2 4 5 1 3]);
                Rs{Si-i} = R;
                Rs{Si-i+1} =[];
                else
                L =ncon({Ls{i-1},MP{i}},{[1,-1,-2],[1,-3,-4]});
                L = ncon({L,conj(MP{i})},{[-1,1,-2,-3],[1,-4,-5]});
                L = ncon({HeisMPO{i},L},{[1,2,-2,3],[1,2,-1,3,-3]});    
                Ls{i} = L;
                R= ncon({MP{i+1},HeisMPO{i+1},conj(MP{i+1})},{[-1,1,-2],[-3,1,-4,2],[-5,2,-6]});
                R = reshape(R,size(R,1),size(R,3),size(R,5));

                Rs{Si-i} = R;
                Rs{Si-i+1} =[];
                    
                end  
        end
     %   x = x-1;
      % [MP,~]=newfDMRG(MP,MO,N,Ls,Rs,x);
        % end for loop
    %elseif(N~=0 && x==0)
    while(N~=0)    
 %The left sweep    
          
          for i= (2*t:-1:2)     
              
              
            if (i==2*t)
                
                
                
                [A1,E]=applyMPOsite(Ls{i-1},HeisMPO{i},[]);
                
                A1 = reshape(A1, [size(A1,1), size(A1,2)*size(A1,3)]);
                [U,S,V] = svd(A1,'econ');
                MP{i} = reshape(V',[size(V,2), size(MP{i},2), size(V,1)/size(MP{i},2)]);
                MP{i-1} = ncon({(MP{i-1}),U*S},{[-1 -2 1],[1 -3]});
                
                L =ncon({Ls{i-2},MP{i-1}},{[1,-1,-2],[1,-3,-4]});
                L = ncon({L,conj(MP{i-1})},{[-1,1,-2,-3],[1,-4,-5]});
                L = ncon({HeisMPO{i-1},L},{[1,2,-2,3],[1,2,-1,3,-3]});    
                Ls{i-1} = L;
                
                
                R= ncon({MP{i},HeisMPO{i},conj(MP{i})},{[-1,1,-2],[-3,1,-4,2],[-5,2,-6]});
                R = reshape(R,size(R,1),size(R,3),size(R,5));
                Rs{1} = R;
                

                
                
                
            elseif(i==2)
                
                [A1,E]=applyMPOsite(Ls{i-1},HeisMPO{i},Rs{Si-i});
                
                A1 = reshape(A1, [size(A1,1), size(A1,2)*size(A1,3)]);
                [U,S,V] = svd(A1,'econ');
                MP{i} = reshape(V',[size(V,2), size(MP{i},2), size(V,1)/size(MP{i},2)]);
                MP{i-1} = ncon({(MP{i-1}),U*S},{[-1 -2 1],[1 -3]});
                
                
       
                
                R= ncon({MP{i},conj(MP{i}),HeisMPO{i},Rs{Si-i}},{[-1,1,2],[-3,3,4],[-2,1,5,3],[2,5,4]},[2 4 5 1 3]);
                

                Rs{Si-i+1} = R;
            else
                [A1,E]=applyMPOsite(Ls{i-1},HeisMPO{i},Rs{Si-i});
                
                A1 = reshape(A1, [size(A1,1), size(A1,2)*size(A1,3)]);
                [U,S,V] = svd(A1,'econ');
                MP{i} = reshape(V',[size(V,2), size(MP{i},2), size(V,1)/size(MP{i},2)]);
                MP{i-1} = ncon({(MP{i-1}),U*S},{[-1 -2 1],[1 -3]});
                
                L =ncon({Ls{i-2},MP{i-1}},{[1,-1,-2],[1,-3,-4]});
                L = ncon({L,conj(MP{i-1})},{[-1,1,-2,-3],[1,-4,-5]});
                L = ncon({HeisMPO{i-1},L},{[1,2,-2,3],[1,2,-1,3,-3]});    
                Ls{i-1} = L;
                
                R= ncon({MP{i},conj(MP{i}),HeisMPO{i},Rs{Si-i}},{[-1,1,2],[-3,3,4],[-2,1,5,3],[2,5,4]},[2 4 5 1 3]);
                

                Rs{Si-i+1} = R;
                
            end             
            
          end
         %The right sweep
          for i= (1:2*t-1)                        
            if (i==1)
                [A1,E]=applyMPOsite([],HeisMPO{1},Rs{Si-1});
                A1 = reshape(A1, [size(A1,1)*size(A1,2), size(A1,3)]); % Convert the rank 3 tensor into a matrix to perform SVD
                [U,S,V] = svd(A1,'econ'); % econ helps in memory 
                MP{1} = reshape(U,[size(U,1)/size(MP{i},2), size(MP{i},2), size(U,2)]);
                MP{i+1} = ncon({(S*V'),MP{i+1}},{[-1 1],[1 -2 -3]});
                
               
                L= ncon({MP{1},HeisMPO{1},conj(MP{1})},{[-1,1,-2],[-3,1,-4,2],[-5,2,-6]});
                L = reshape(L,size(L,6),size(L,4),size(L,2));
                Ls{1} =L;
                
                
                R= ncon({MP{i+1},conj(MP{i+1}),HeisMPO{i+1},Rs{Si-i-1}},{[-1,1,2],[-3,3,4],[-2,1,5,3],[2,5,4]},[2 4 5 1 3]);
                Rs{Si-i} = R;
                Rs{Si-i+1} =[];
                
                
                
            elseif (i<2*t-1)
                [A1,E]=applyMPOsite(Ls{i-1},HeisMPO{i},Rs{Si-i});
                A1 = reshape(A1, [size(A1,1)*size(A1,2), size(A1,3)]); % Convert the rank 3 tensor into a matrix to perform SVD
                [U,S,V] = svd(A1,'econ'); % econ helps in memory 
                MP{i} = reshape(U,[size(U,1)/size(MP{i},2), size(MP{i},2), size(U,2)]);
                MP{i+1} = ncon({(S*V'),MP{i+1}},{[-1 1],[1 -2 -3]});
                
                L =ncon({Ls{i-1},MP{i}},{[1,-1,-2],[1,-3,-4]});
                L = ncon({L,conj(MP{i})},{[-1,1,-2,-3],[1,-4,-5]});
                L = ncon({HeisMPO{i},L},{[1,2,-2,3],[1,2,-1,3,-3]});    
                Ls{i} = L;
                
                R= ncon({MP{i+1},conj(MP{i+1}),HeisMPO{i+1},Rs{Si-i-1}},{[-1,1,2],[-3,3,4],[-2,1,5,3],[2,5,4]},[2 4 5 1 3]);
                Rs{Si-i} = R;
                Rs{Si-i+1} =[];
            else
                [A1,E]=applyMPOsite(Ls{i-1},HeisMPO{i},Rs{Si-i});
                A1 = reshape(A1, [size(A1,1)*size(A1,2), size(A1,3)]); % Convert the rank 3 tensor into a matrix to perform SVD
                [U,S,V] = svd(A1,'econ'); % econ helps in memory 
                MP{i} = reshape(U,[size(U,1)/size(MP{i},2), size(MP{i},2), size(U,2)]);
                MP{i+1} = ncon({(S*V'),MP{i+1}},{[-1 1],[1 -2 -3]});
                
                L =ncon({Ls{i-1},MP{i}},{[1,-1,-2],[1,-3,-4]});
                L = ncon({L,conj(MP{i})},{[-1,1,-2,-3],[1,-4,-5]});
                L = ncon({HeisMPO{i},L},{[1,2,-2,3],[1,2,-1,3,-3]});    
                Ls{i} = L;
                
                R= ncon({MP{i+1},HeisMPO{i+1},conj(MP{i+1})},{[-1,1,-2],[-3,1,-4,2],[-5,2,-6]});
                R = reshape(R,size(R,1),size(R,3),size(R,5));

                Rs{Si-i} = R;
                Rs{Si-i+1} =[];
                    
                                                    
            end             
            
         end
          N=N-1;
          
         %[MP,~]=newfDMRG(MP,MO,N,Ls,Rs,x);
         
    
    
        
       

    
    
    
    end

 A=MP;
 %E= calcR(MP(1:2*t),HeisMPO(1:2*t));
 
end