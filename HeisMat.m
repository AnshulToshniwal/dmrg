function [A] = HeisMat(s) 
global J1;
[Sz,Sp,Sm,I] = Spin(s);
A= J1*((kron(Sp,I)*kron(I,Sm)+kron(Sm,I)*kron(I,Sp))/2 + (kron(Sz,I)*kron(I,Sz)));
end
