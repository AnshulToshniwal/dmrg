function [A] = MPOP(N,J2,s)

[Sz,Sp,Sm,I] = Spin(s); % get local spin operators

Hloc = cell(11,11);
Hloc(:) = {zeros(size(I))};
Hloc{1,1} = I;
Hloc{2,1} = (Sm*Sm);
Hloc{3,1} = (Sp*Sp);
Hloc{4,1} = (Sz*Sz);
Hloc{5,1} = (Sm*Sp);
Hloc{6,1} = (Sm*Sz);
Hloc{7,1} = (Sp*Sm);
Hloc{8,1} = (Sp*Sz);
Hloc{9,1} = (Sz*Sp);
Hloc{10,1} = (Sz*Sm);

Hloc{11,2} = (J2/4)*(Sp*Sp);
Hloc{11,3} = (J2/4)*(Sm*Sm);
Hloc{11,4} = (J2)*(Sz*Sz);
Hloc{11,5} = (J2/4)*(Sp*Sm);
Hloc{11,6} = (J2/2)*(Sp*Sz);
Hloc{11,7} = (J2/4)*(Sm*Sp);
Hloc{11,8} = (J2/2)*(Sm*Sz);
Hloc{11,9} = (J2/2)*(Sz*Sm);
Hloc{11,10}= (J2/2)*(Sz*Sp);
Hloc{11,11}= I;


%reshape the cell into an array for tensor of rank 4 
Hloc = cell2mat(reshape(Hloc,[1 1 size(Hloc,1) size(Hloc,2)])); 
Hloc = permute(Hloc,[3 2 4 1]); % leg order: left-bottom-right-top


% full chain
Hs = cell(1,N); %Total hamiltonian
Hs(:) = {Hloc};%Initialising
Hs{1} = Hs{1}(11,:,:,:); % choose the last components of the left leg as first tensor is rank 3 tensor
Hs{N} = Hs{N}(:,:,1,:);
A=Hs;
end