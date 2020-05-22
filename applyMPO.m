function [A1,A2,E] = applyMPO(Hl,Hloc1,Hloc2,Hr,s)
 
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
global D;
Hn1 = ncon({Hloc1,Hloc2},{[-1,-2,1,-3],[1,-4,-5,-6]}); % checked !
Hn2 = ncon({Hl,Hn1,Hr},{[-1,1,-2],[1,-3,-4,-5,2,-6],[-7,2,-8]}); %checked twice of the order

Hn2t = permute(Hn2,[1 3 5 7 2 4 6 8]);
t = size(Hn2t);
Hn3 = reshape(Hn2t,[(t(1)*t(2)*t(3)*t(4)),(t(5)*t(6)*t(7)*t(8))]);
r = rand(size(Hn3,1),1);

[V1,E1] = modifiedlanczos(Hn3,r); % hope it doesnt make a mess

%[V1,E1] = eigs(Hn3,1,'sr'); 
E= E1;
%need to reshape the matrix accordin to local space dimension
alpha = size(Hl,1);
beta = size(Hr,1);

%t =sqrt(size(V1,1));
V = reshape(V1,[alpha*(2*s+1),beta*(2*s+1)]);
[U2,S2,V2] = svd(V);
X = V2';
if ((size(S2,1))>D)
    U2 = U2(:,1:D);
    S2= S2(1:D,1:D);
    X = X(1:D,:);
    A1 = reshape(U2,alpha,(2*s+1),[]);
    A2 = reshape(X,[],(2*s+1),beta);
else
    A1 = reshape(U2,alpha,(2*s+1),[]);
    A2 = reshape(X,[],(2*s+1),beta);
    
      
end

end
