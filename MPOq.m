function [A] = MPOq(N,J1,J2,s)


[Sz,Sp,Sm,I] = Spin(s); % get local spin operators

Hloc = cell(14,14);

Hloc(:) = {zeros(size(I))};
Hloc{1,1} = I;
Hloc{2,1} = (Sm^2);
Hloc{3,1} = (Sp^2);
Hloc{4,1} = (Sz^2);
Hloc{5,1} = (Sm*Sp);
Hloc{6,1} = (Sm*Sz);
Hloc{7,1} = (Sp*Sm);
Hloc{8,1} = (Sp*Sz);
Hloc{9,1} = (Sz*Sp);
Hloc{10,1} = (Sz*Sm);
Hloc{11,1} = (Sm);
Hloc{12,1} = (Sp);
Hloc{13,1} = (Sz);

Hloc{14,2} = (J2/4)*(Sp^2);
Hloc{14,3} = (J2/4)*(Sm^2);
Hloc{14,4} = J2*(Sz^2);
Hloc{14,5} = (J2/4)*(Sp*Sm);
Hloc{14,6} = (J2/2)*(Sp*Sz);
Hloc{14,7} = (J2/4)*(Sm*Sp);
Hloc{14,8} = (J2/2)*(Sm*Sz);
Hloc{14,9} = (J2/2)*(Sz*Sm);
Hloc{14,10}= (J2/2)*(Sz*Sp);
Hloc{14,11}= (J1/2)*(Sp);
Hloc{14,12}= (J1/2)*(Sm);
Hloc{14,13}= (J1)*(Sz);
Hloc{14,14}= I;

%reshape the cell into an array for tensor of rank 4 
Hloc = cell2mat(reshape(Hloc,[1 1 size(Hloc,1) size(Hloc,2)])); 
Hloc = permute(Hloc,[3 2 4 1]); % leg order: left-bottom-right-top


% full chain
Hs = cell(1,N); %Total hamiltonian
Hs(:) = {Hloc};%Initialising
Hs{1} = Hs{1}(14,:,:,:); % choose the last components of the left leg as first tensor is rank 3 tensor
Hs{N} = Hs{N}(:,:,1,:);
A=Hs;
end