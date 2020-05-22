function [A1,E] = applyMPO1site(Hl,Hloc,Hr,s)
 
%Hl, Hloc1,Hloc2, Hr: [tensors] The Hamiltonian for the left, 2 sites , and
%       right parts of the chain.Hloc1 and hloc2 do not depend on MPS but
%       Hl and Hr do.Hl and Hr are rank 3 tensors while Hloc1 and Hloc2 are
%       rank 4 tensors
% Ai : [tensor] Ket tensor as an initial vector for the Lanczos method
%< Output >
% A1,A2 : [tensor] Ket tensors as the ground state of the effective
%       Hamiltonian.
% E : [numeric] Ground-state energy of the  Hamiltonian.
% s: spin of local space
if isempty(Hl) && isempty(Hr)
    error('ERR: Hleft and Hright cannot be empty at the same time.');
end
if isempty(Hl)
    
    Hloc = reshape(Hloc,[1 size(Hloc)]);
    Hn = ncon({Hloc,Hr},{[-1,-2,-3,1,-4],[-5,1,-6]});
    Hn2 = permute(Hn,[1 3 5 2 4 6]);
    t = size(Hn2);
    Hn3 = reshape(Hn2,[(t(1)*t(2)*t(3)),(t(4)*t(5)*t(6))]);

 
   [V1,E1] = eigs(Hn3,1,'sr'); 
    E= E1;
    A1 = reshape(V1,1,(2*s+1),[]);
elseif isempty(Hr) %need to check
    Hn2 = ncon({Hl,Hloc},{[-1,1,-2],[1,-3,-4,-5]});
    Hn = permute(Hn2,[1 2 3 5 4]);
    Hn2 = permute(Hn,[1 3 5 2 4 6]);
    t = size(Hn2);
    Hn3 = reshape(Hn2,[(t(1)*t(3)*t(5)),(t(2)*t(4)*t(6))]);


    [V1,E1] = eigs(Hn3,1,'sr'); 
     E= E1;
     A1 = reshape(V1,[],(2*s+1),1);
    
else

global D;
Hn2 = ncon({Hl,Hloc,Hr},{[-1,1,-2],[1,-3,2,-4],[-5,2,-6]}); %checked twice of the order

Hn2t = permute(Hn2,[1 3 5 2 4 6]);
t = size(Hn2t);
Hn3 = reshape(Hn2t,[(t(1)*t(2)*t(3)),(t(4)*t(5)*t(6))]);

r = rand(size(Hn3,1),1);

%[V1,E1,~] = lanczos(Hn3,r); % hope it doesnt make a mess

[V1,E1] = eigs(Hn3,1,'sr'); 
E= E1;
%need to reshape the matrix accordin to local space dimension
alpha = size(Hl,1);
beta = size(Hr,1);


A1 = reshape(V1,[alpha,(2*s+1),beta]);


end
end