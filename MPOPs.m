function [A] = MPOPs(N,J1,J2,J3,s)

[Sz,Sp,Sm,I] = Spin(s); % get local spin operators


Hs1 = cell(45,45);
Hs1(:) = {zeros(size(I))};
Hs1{1,1} = I;
Hs1{2,1} = (Sp);
Hs1{3,1} = (Sm);
Hs1{4,1} = (Sz);
Hs1{5,2} = (J1/2)*(Sm);
Hs1{5,3} = (J1/2)*(Sp);
Hs1{5,4}= J1*Sz;
Hs1{5,5} = I;
Hs1{6,6} = I;
Hs1{7,6} = (Sm^2);
Hs1{8,6} = (Sp^2);
Hs1{9,6} = (Sz^2);
Hs1{10,6} = (Sm*Sp);
Hs1{11,6} = (Sm*Sz);
Hs1{12,6} = (Sp*Sm);
Hs1{13,6} = (Sp*Sz);
Hs1{14,6} = (Sz*Sp);
Hs1{15,6} = (Sz*Sm);

Hs1{16,7} = (J2/4)*(Sp^2);
Hs1{16,8} = (J2/4)*(Sm^2);
Hs1{16,9} = J2*(Sz^2);
Hs1{16,10} = (J2/4)*(Sp*Sm);
Hs1{16,11} = (J2/2)*(Sp*Sz);
Hs1{16,12} = (J2/2)*(Sm*Sp);
Hs1{16,13} = (J2/2)*(Sm*Sz);
Hs1{16,14} = (J2/2)*(Sz*Sm);
Hs1{16,15}= (J2/2)*(Sz*Sp);
Hs1{16,16}= I;

Hs1{17,17} = I;
Hs1{18,17} = (Sm^3);
Hs1{19,17} = ((Sm^2)*Sp);
Hs1{20,17} = ((Sm^2)*Sz);
Hs1{21,17} = (Sm*Sz*Sm);
Hs1{22,17} = (Sm*Sz*Sp);
Hs1{23,17} = (Sm*(Sz^2));
Hs1{24,17} = (Sp*(Sm^2));
Hs1{25,17} = (Sp*(Sm*Sp));
Hs1{26,17} = (Sp*Sm*Sz);
Hs1{27,17} = ((Sp^2)*Sm);
Hs1{28,17} = (Sp^3);
Hs1{29,17} = ((Sp^2)*Sz);
Hs1{30,17} = ((Sp)*Sz*(Sm));
Hs1{31,17} = (Sp*Sz*Sp);
Hs1{32,17} = (Sp*(Sz^2));
Hs1{33,17} = (Sz*(Sm^2));
Hs1{34,17} = (Sz*(Sm*Sp));
Hs1{35,17} = (Sz*(Sm*Sz));
Hs1{36,17} = (Sz*Sp*Sm);
Hs1{37,17} = (Sz*(Sp^2));
Hs1{38,17} = (Sz*Sp*Sz);
Hs1{39,17} = ((Sz^2)*Sm);
Hs1{40,17} = ((Sz^2)*Sp);
Hs1{41,17} = (Sz^3);
Hs1{42,17} = (Sm*Sp*Sm);
Hs1{43,17} = (Sm*(Sp^2));
Hs1{44,17} = (Sm*Sp*Sz);

Hs1{45,18} = (J3/8)*(Sp^3);
Hs1{45,19} = (J3/8)*((Sp^2)*Sm);
Hs1{45,20} = (J3/4)*(Sp^2)*Sz;
Hs1{45,21} = (J3/4)*(Sp*Sz*Sp);
Hs1{45,22} = (J3/4)*(Sp*Sz*Sm);
Hs1{45,23} = (J3/2)*(Sp*(Sz^2));
Hs1{45,24} = (J3/8)*(Sm*(Sp^2));
Hs1{45,25} = (J3/8)*(Sm*Sp*Sm);
Hs1{45,26}= (J3/4)*(Sm*Sp*Sz);
Hs1{45,27}= (J3/8)*((Sm^2)*(Sp));
Hs1{45,28} = (J3/8)*(Sm^3);
Hs1{45,29} = (J3/4)*(Sm^2)*Sz;
Hs1{45,30} = (J3/4)*(Sm*Sz*Sp);
Hs1{45,31} = (J3/4)*(Sm*Sz*Sm);
Hs1{45,32} = (J3/2)*(Sm*(Sz^2));
Hs1{45,33} = (J3/4)*(Sz*(Sp^2));
Hs1{45,34} = (J3/4)*(Sz*Sp*Sm);
Hs1{45,35} = (J3/2)*(Sz*Sp*Sz);
Hs1{45,36}= (J3/4)*(Sz*Sm*Sp);
Hs1{45,37}= (J3/4)*(Sz*(Sm^2));
Hs1{45,38} = (J3/2)*(Sz*Sm*Sz);
Hs1{45,39} = (J3/2)*((Sz^2)*Sp);
Hs1{45,40} = (J3/2)*((Sz^2)*Sm);
Hs1{45,41} = (J3)*(Sz^3);
Hs1{45,42}= (J3/8)*(Sp*Sm*Sp);
Hs1{45,43} = (J3/8)*(Sp*(Sm^2));
Hs1{45,44} = (J3/4)*(Sp*Sm*Sz);
Hs1{45,45} = I;
%Setting the first and last values

Hf = cell(1,45);
Hf(:) = {zeros(size(I))};

Hf{1,2} = (J1/2)*(Sm);
Hf{1,3} = (J1/2)*(Sp);
Hf{1,4}= J1*Sz;
Hf{1,5} = I;
Hf{1,7} = (J2/4)*(Sp^2);
Hf{1,8} = (J2/4)*(Sm^2);
Hf{1,9} = J2*(Sz^2);
Hf{1,10} = (J2/4)*(Sp*Sm);
Hf{1,11} = (J2/2)*(Sp*Sz);
Hf{1,12} = (J2/2)*(Sm*Sp);
Hf{1,13} = (J2/2)*(Sm*Sz);
Hf{1,14} = (J2/2)*(Sz*Sm);
Hf{1,15}= (J2/2)*(Sz*Sp);
Hf{1,16}= I;
Hf{1,18} = (J3/8)*(Sp^3);
Hf{1,19} = (J3/8)*((Sp^2)*Sm);
Hf{1,20} = (J3/4)*(Sp^2)*Sz;
Hf{1,21} = (J3/4)*(Sp*Sz*Sp);
Hf{1,22} = (J3/4)*(Sp*Sz*Sm);
Hf{1,23} = (J3/2)*(Sp*(Sz^2));
Hf{1,24} = (J3/8)*(Sm*(Sp^2));
Hf{1,25} = (J3/8)*(Sm*Sp*Sm);
Hf{1,26}= (J3/4)*(Sm*Sp*Sz);
Hf{1,27}= (J3/8)*((Sm^2)*(Sp));
Hf{1,28} = (J3/8)*(Sm^3);
Hf{1,29} = (J3/4)*(Sm^2)*Sz;
Hf{1,30} = (J3/4)*(Sm*Sz*Sp);
Hf{1,31} = (J3/4)*(Sm*Sz*Sm);
Hf{1,32} = (J3/2)*(Sm*(Sz^2));
Hf{1,33} = (J3/4)*(Sz*(Sp^2));
Hf{1,34} = (J3/4)*(Sz*Sp*Sm);
Hf{1,35} = (J3/2)*(Sz*Sp*Sz);
Hf{1,36}= (J3/4)*(Sz*Sm*Sp);
Hf{1,37}= (J3/4)*(Sz*(Sm^2));
Hf{1,38} = (J3/2)*(Sz*Sm*Sz);
Hf{1,39} = (J3/2)*((Sz^2)*Sp);
Hf{1,40} = (J3/2)*((Sz^2)*Sm);
Hf{1,41} = (J3)*(Sz^3);
Hf{1,42}= (J3/8)*(Sp*Sm*Sp);
Hf{1,43} = (J3/8)*(Sp*(Sm^2));
Hf{1,44} = (J3/4)*(Sp*Sm*Sz);
Hf{1,45} = I;


Hl = cell(45,1);
Hl(:) = {zeros(size(I))};
Hl{1,1} = I;
Hl{2,1} = (Sp);
Hl{3,1} = (Sm);
Hl{4,1} = (Sz);

Hl{6,1} = I;
Hl{7,1} = (Sm^2);
Hl{8,1} = (Sp^2);
Hl{9,1} = (Sz^2);
Hl{10,1} = (Sm*Sp);
Hl{11,1} = (Sm*Sz);
Hl{12,1} = (Sp*Sm);
Hl{13,1} = (Sp*Sz);
Hl{14,1} = (Sz*Sp);
Hl{15,1} = (Sz*Sm);

Hl{17,1} = I;
Hl{18,1} = (Sm^3);
Hl{19,1} = ((Sm^2)*Sp);
Hl{20,1} = ((Sm^2)*Sz);
Hl{21,1} = (Sm*Sz*Sm);
Hl{22,1} = (Sm*Sz*Sp);
Hl{23,1} = (Sm*(Sz^2));
Hl{24,1} = (Sp*(Sm^2));
Hl{25,1} = (Sp*(Sm*Sp));
Hl{26,1} = (Sp*Sm*Sz);
Hl{27,1} = ((Sp^2)*Sm);
Hl{28,1} = (Sp^3);
Hl{29,1} = ((Sp^2)*Sz);
Hl{30,1} = ((Sp)*Sz*(Sm));
Hl{31,1} = (Sp*Sz*Sp);
Hl{32,1} = (Sp*(Sz^2));
Hl{33,1} = (Sz*(Sm^2));
Hl{34,1} = (Sz*(Sm*Sp));
Hl{35,1} = (Sz*(Sm*Sz));
Hl{36,1} = (Sz*Sp*Sm);
Hl{37,1} = (Sz*(Sp^2));
Hl{38,1} = (Sz*Sp*Sz);
Hl{39,1} = ((Sz^2)*Sm);
Hl{40,1} = ((Sz^2)*Sp);
Hl{41,1} = (Sz^3);
Hl{42,1} = (Sm*Sp*Sm);
Hl{43,1} = (Sm*(Sp^2));
Hl{44,1} = (Sm*Sp*Sz);

%reshape the cell into an array for tensor of rank 4 
Hf1 = cell2mat(reshape(Hf,[1 1 size(Hf,1) size(Hf,2)])); 
Hf1 = permute(Hf1,[3 2 4 1]); % leg order: left-bottom-right-top

Hl1 = cell2mat(reshape(Hl,[1 1 size(Hl,1) size(Hl,2)])); 
Hl1 = permute(Hl1,[3 2 4 1]); % leg order: left-bottom-right-top

Hs1 = cell2mat(reshape(Hs1,[1 1 size(Hs1,1) size(Hs1,2)])); 
Hs1 = permute(Hs1,[3 2 4 1]); % leg order: left-bottom-right-top






% full chain
Hs = cell(1,N); %Total hamiltonian
Hs(:) = {Hs1};%Initialising
Hs(1) = {Hf1}; % choose the last components of the left leg as first tensor is rank 3 tensor
Hs(N) = {Hl1};
A=Hs;
end