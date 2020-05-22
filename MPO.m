% construct MPO for heisenberg hamiltonian
function [A] = MPO(N,J1,s)

[Sz,Sp,Sm,I] = Spin(s); % get local spin operators

Hloc = cell(5,5);
Hloc(:) = {zeros(size(I))};
Hloc{1,1} = I;
Hloc{2,1} = (Sp);
Hloc{3,1} = (Sm);
Hloc{4,1} = (Sz);
Hloc{5,2} = (J1/2)*(Sm);
Hloc{5,3} = (J1/2)*(Sp);
Hloc{5,4}= J1*Sz;
Hloc{5,5} = I;

%reshape the cell into an array for tensor of rank 4 
Hloc = cell2mat(reshape(Hloc,[1 1 size(Hloc,1) size(Hloc,2)])); 
Hloc = permute(Hloc,[3 2 4 1]); % leg order: left-bottom-right-top


% full chain
Hs = cell(1,N); %Total hamiltonian
Hs(:) = {Hloc};%Initialising
Hs{1} = Hs{1}(5,:,:,:); % choose the last components of the left leg as first tensor is rank 3 tensor
Hs{N} = Hs{N}(:,:,1,:);
A=Hs;
end