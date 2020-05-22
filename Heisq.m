function [A] = Heisq(s) 
global J1 J2 ;
[Sz,Sp,Sm,I] = Spin(s);
A= J1*((kron(Sp,I)*kron(I,Sm)+kron(Sm,I)*kron(I,Sp))/2 + (kron(Sz,I)*kron(I,Sz)))....
 + J2*(((kron(Sp,I)*kron(I,Sm)+kron(Sm,I)*kron(I,Sp))/2 + (kron(Sz,I)*kron(I,Sz)))^2);
end