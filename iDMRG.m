function [A,E,x,y] = iDMRG(Hamil,N,s)
%Hl, Hloc1,Hloc2, Hr: [tensors] The Hamiltonian for the left, 2 sites , and
%       right parts of the chain.Hloc1 and hloc2 do not depend on MPS but
%       Hl and Hr do.Hl and Hr are rank 3 tensors while Hloc1 and Hloc2 are
%       rank 4 tensors
% Ai : an initial vector for the Lanczos method
% s: spin of local space
% N: 2*number of sites
%< Output >
% A1 : [tensor] MPS of the ground state
% E : [numeric] Ground-state energy of the  Hamiltonian.
%
global J1 J2 J3

Si = 2*N; % Number of sites 
M = cell(1,Si); %the ground state MPS 
Ls= cell(1,Si);
Rs= cell(1,Si);
if (strcmp(Hamil,'Heisenberg'))
    
    % Performing two site calculation exactly for spin 1/2
    %Heis = [J1/4 0 0 0; 0 -J1/4 J1/2 0; 0 J1/2 -J1/4 0; 0 0 0 J1/4];
    Heis = HeisMat(s);
    [V,~] = eigs(Heis,1,'sa');
    V =V/norm(V);

     % reshaping according to local space dimension
    dim = 2*s+1; 
    V = reshape(V,dim,dim);
    [U1,S1,V1] = svd(V);
    V1 = V1';
    U1 = reshape(U1,1,dim,[]);
    V1 = reshape(V1,[],dim,1);
    M{1} = U1;
    M{Si} = V1;
    
    HeisMPOt = MPO(Si,J1,s);
    
    L = calcL(M(1:1),HeisMPOt(1:1));
    R= calcR(M(Si:Si),HeisMPOt(Si:Si));
    Ls{1} = L;
    Rs{1} = R;
    for i= (2:N)
        %{} in cell returns a tensor while () returns cell
        %HeisMPO = MPO(2*i,J1,s);
        if (i==2)
            [A1,A2,E1]=applyMPO(L,HeisMPOt{i},HeisMPOt{i+1},R,s);
            M{i} = A1;
            M{Si-i+1} = A2;
        else
        L =ncon({L,M{i-1}},{[1,-1,-2],[1,-3,-4]});
        L = ncon({L,conj(M{i-1})},{[-1,1,-2,-3],[1,-4,-5]});
        L = ncon({HeisMPOt{i-1},L},{[1,2,-2,3],[1,2,-1,3,-3]});    
        Ls{i-1} = L;
        R= ncon({M{Si-i+2},conj(M{Si-i+2}),HeisMPOt{Si-i+2},R},{[-1,1,2],[-3,3,4],[-2,1,5,3],[2,5,4]},[2 4 5 1 3]);
        Rs{i-1} = R;
        [A1,A2,E1]=applyMPO(L,HeisMPOt{i},HeisMPOt{i+1},R,s);
        M{i} = A1;
        M{Si-i+1} = A2;
        end
        
        
    end
    Last =ncon({Ls{N-1},M{N}},{[1,-1,-2],[1,-3,-4]});
    Last = ncon({Last,conj(M{N})},{[-1,1,-2,-3],[1,-4,-5]});
    Last = ncon({HeisMPOt{N},Last},{[1,2,-2,3],[1,2,-1,3,-3]});    
    Ls{N} = Last;
    Rst= ncon({M{N+1},conj(M{N+1}),HeisMPOt{N+1},R},{[-1,1,2],[-3,3,4],[-2,1,5,3],[2,5,4]},[2 4 5 1 3]);
    Rs{N} = Rst;
    
    A = M;
    E=E1;
    x = Ls;
    y= Rs;
end
if (strcmp(Hamil,'Heisenberg_sumpower'))
    
   Heisp = HeisP(s);
    [V,~] = eigs(Heisp,1,'sr');
    V=V/norm(V);
    % reshaping according to local space dimension
    dim = 2*s+1; 
    V = reshape(V,dim,dim);
    [U1,S1,V1] = svd(V);
    V1 = V1';
    U1 = reshape(U1,1,dim,[]);
    V1 = reshape(V1,[],dim,1);
    
    M{1} = U1;
    M{Si} = V1;
    
    HeisMPOt = test(Si,J1,J2,J3,s);
    
    L = calcL(M(1:1),HeisMPOt(1:1));
    R= calcR(M(Si:Si),HeisMPOt(Si:Si));
    Ls{1} =L;
    Rs{1} = R;
    for i= (2:N)
        %{} in cell returns a tensor while () returns cell
        
        if (i==2)
            [A1,A2,E1]=applyMPO(L,HeisMPOt{i},HeisMPOt{i+1},R,s);
            M{i} = A1;
            M{Si-i+1} = A2;
        else
        L =ncon({L,M{i-1}},{[1,-1,-2],[1,-3,-4]});
        L = ncon({L,conj(M{i-1})},{[-1,1,-2,-3],[1,-4,-5]});
        L = ncon({HeisMPOt{i-1},L},{[1,2,-2,3],[1,2,-1,3,-3]}); 
        Ls{i-1} = L;
        R= ncon({M{Si-i+2},conj(M{Si-i+2}),HeisMPOt{Si-i+2},R},{[-1,1,2],[-3,3,4],[-2,1,5,3],[2,5,4]},[2 4 5 1 3]);
        Rs{i-1} = R;
        [A1,A2,E1]=applyMPO(L,HeisMPOt{i},HeisMPOt{i+1},R,s);
        M{i} = A1;
        M{Si-i+1} = A2;
        end
        
    end
    Last =ncon({L,M{N}},{[1,-1,-2],[1,-3,-4]});
    Last = ncon({Last,conj(M{N})},{[-1,1,-2,-3],[1,-4,-5]});
    Last = ncon({HeisMPOt{N},Last},{[1,2,-2,3],[1,2,-1,3,-3]});    
    Ls{N} = Last;
    Rst= ncon({M{N+1},conj(M{N+1}),HeisMPOt{N+1},R},{[-1,1,2],[-3,3,4],[-2,1,5,3],[2,5,4]},[2 4 5 1 3]);
    Rs{N} = Rst;
    
    A = M;
    E=E1;
    x= Ls;
    y= Rs;
end
end