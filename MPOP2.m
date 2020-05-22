function [A] = MPOP2(N,J3,s)

[Sz,Sp,Sm,I] = Spin(s); % get local spin operators

Hloc = cell(29,29);
Hloc(:) = {zeros(size(I))};
Hloc{1,1} = I;
Hloc{2,1} = (Sm^3);
Hloc{3,1} = ((Sm^2)*Sp);
Hloc{4,1} = ((Sm^2)*Sz);
Hloc{5,1} = (Sm*Sz*Sm);
Hloc{6,1} = (Sm*Sz*Sp);
Hloc{7,1} = (Sm*(Sz^2));
Hloc{8,1} = (Sp*(Sm^2));
Hloc{9,1} = (Sp*(Sm*Sp));
Hloc{10,1} = (Sp*Sm*Sz);
Hloc{11,1} = ((Sp^2)*Sm);
Hloc{12,1} = (Sp^3);
Hloc{13,1} = ((Sp^2)*Sz);
Hloc{14,1} = ((Sp)*Sz*(Sm));
Hloc{15,1} = (Sp*Sz*Sp);
Hloc{16,1} = (Sp*(Sz^2));
Hloc{17,1} = (Sz*(Sm^2));
Hloc{18,1} = (Sz*(Sm*Sp));
Hloc{19,1} = (Sz*(Sm*Sz));
Hloc{20,1} = (Sz*Sp*Sm);
Hloc{21,1} = (Sz*(Sp^2));
Hloc{22,1} = (Sz*Sp*Sz);
Hloc{23,1} = ((Sz^2)*Sm);
Hloc{24,1} = ((Sz^2)*Sp);
Hloc{25,1} = (Sz^3);
Hloc{26,1} = (Sm*Sp*Sm);
Hloc{27,1} = (Sm*(Sp^2));
Hloc{28,1} = (Sm*Sp*Sz);

Hloc{29,2} = (J3/8)*(Sp^3);
Hloc{29,3} = (J3/8)*((Sp^2)*Sm);
Hloc{29,4} = (J3/4)*(Sp^2)*Sz;
Hloc{29,5} = (J3/4)*(Sp*Sz*Sp);
Hloc{29,6} = (J3/4)*(Sp*Sz*Sm);
Hloc{29,7} = (J3/2)*(Sp*(Sz^2));
Hloc{29,8} = (J3/8)*(Sm*(Sp^2));
Hloc{29,9} = (J3/8)*(Sm*Sp*Sm);
Hloc{29,10}= (J3/4)*(Sm*Sp*Sz);
Hloc{29,11}= (J3/8)*((Sm^2)*(Sp));
Hloc{29,12} = (J3/8)*(Sm^3);
Hloc{29,13} = (J3/4)*(Sm^2)*Sz;
Hloc{29,14} = (J3/4)*(Sm*Sz*Sp);
Hloc{29,15} = (J3/4)*(Sm*Sz*Sm);
Hloc{29,16} = (J3/2)*(Sm*(Sz^2));
Hloc{29,17} = (J3/4)*(Sz*(Sp^2));
Hloc{29,18} = (J3/4)*(Sz*Sp*Sm);
Hloc{29,19} = (J3/2)*(Sz*Sp*Sz);
Hloc{29,20}= (J3/4)*(Sz*Sm*Sp);
Hloc{29,21}= (J3/4)*(Sz*(Sm^2));
Hloc{29,22} = (J3/2)*(Sz*Sm*Sz);
Hloc{29,23} = (J3/2)*((Sz^2)*Sp);
Hloc{29,24} = (J3/2)*((Sz^2)*Sm);
Hloc{29,25} = (J3)*(Sz^3);
Hloc{29,26}= (J3/8)*(Sp*Sm*Sp);
Hloc{29,27} = (J3/8)*(Sp*(Sm^2));
Hloc{29,28} = (J3/4)*(Sp*Sm*Sz);
Hloc{29,29} = I;


%reshape the cell into an array for tensor of rank 4 
Hloc = cell2mat(reshape(Hloc,[1 1 size(Hloc,1) size(Hloc,2)])); 
Hloc = permute(Hloc,[3 2 4 1]); % leg order: left-bottom-right-top


% full chain
Hs = cell(1,N); %Total hamiltonian
Hs(:) = {Hloc};%Initialising
Hs{1} = Hs{1}(29,:,:,:); % choose the last components of the left leg as first tensor is rank 3 tensor
Hs{N} = Hs{N}(:,:,1,:);
A=Hs;
end