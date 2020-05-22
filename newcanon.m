function [A] = newcanon(M,j) %converts MPS in site canonical form
%input check
if ((size(M,2) < j))
    error('ERROR: the site specified is beyond the dimension of MPS.');
end
N = size(M,2);
% end points
if ( j ==1)
    A = rightcanon(M);
elseif(j==N)
    A=leftcanon(M);
else
    A1 = leftcanon(M(1:j));
    A2= M(j:N);
    A2{1}=A1{j};
    A2 = rightcanon(A2);
    A(1:j-1) = A1(1:j-1);
    A(j:N) = A2(1:(N-j+1));
    
end
