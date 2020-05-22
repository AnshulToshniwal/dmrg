function [A] = HeisC(s) 
global  J3 ;
[Sz,Sp,Sm,I] = Spin(s);
A= J3*(((kron(Sp,I)*kron(I,Sm)+kron(Sm,I)*kron(I,Sp))/2 + (kron(Sz,I)*kron(I,Sz)))^3);
end