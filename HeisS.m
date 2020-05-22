function [A] = HeisS(s) 
global  J2 ;
[Sz,Sp,Sm,~] = Spin(s);
A= J2*(((kron(Sp,Sm)+kron(Sm,Sp))/2 + (kron(Sz,Sz)))^2);
end