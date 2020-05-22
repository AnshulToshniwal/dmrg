% Defining generic spin matrices for a given s

function [Sz,Sp,Sm,I] = Spin(s)
l = 2*s+1;
%Sx = zeros(l,l);
%Sy = zeros(l,l);
I = eye(l,l);
Sz = zeros(l,l);
Sp = zeros(l,l);
Sm = zeros(l,l);
for i=1:l
    
    
    for j =1:l
        
        
 %       Sx(i,j) = (kroneckerDelta(sym(i-1-s),sym(j-s))+kroneckerDelta(sym(i-s),sym(j-1-s)))*sqrt(s*(s+1)-(i-1-s)*(j-1-s))/2;
  %      Sy(i,j) = ((kroneckerDelta(sym(i-1-s),sym(j-s))-kroneckerDelta(sym(i-s),sym(j-1-s)))*sqrt(s*(s+1)-(i-s-1)*(j-s-1)))/(2*complex(0,-1));
        Sz(i,j) = -kroneckerDelta(sym(i-1-s),sym(j-1-s))*(j-1-s);
        Sm(i,j) = kroneckerDelta(sym(i-1-s),sym(j-s))*sqrt(s*(s+1)-(i-1-s)*(j-1-s));
        Sp(i,j) = kroneckerDelta(sym(i-s),sym(j-1-s))*sqrt(s*(s+1)-(i-1-s)*(j-1-s));
    end 
end
    