%Lanczos function taking as input the matrix(A), initial seed vector(r),number of iterations(m)
%Outputs the tridiagonal matrix(t), the orthonormal vectors (V) and the
%betas and the alphas
function [V,E,T] = Lanczos(A,r)
if ((size(A,1) ~= size(A,2)))
    error('ERROR: A is not square matrix.');
elseif size(A,1) ~= numel(r)
    error('ERROR: r is not compatible with A.');
end
% m is the number of lanzos vectors considered for the diagonalisation

A= (A+A')/2; % To reduce numerical noise
global Ld;
m=Ld;
alph = zeros(1,m);% initialising the vector
T = zeros(m,m);% initialising the vector
bet = zeros(1,m);% initialising the vector
[s1,~] = size(A);
V = zeros(s1,m);% initialising the vectors
r= r/norm(r);% normalising the seed vector
W = zeros(s1,m);  % initialising the matrix
alph(1) = (r')*(A')*r; % Setting the first value
V(:,1) = r;% Setting the first value
W(:,1) = A*r - alph(1)*r;% Setting the first value

if m<=s1
    for i= 2:m
    bet(i) = norm(W(:,i-1));
    V(:,i) = (W(:,i-1))/ bet(i); 
    alph(i) = (V(:,i)')*(A')*(V(:,i));
    W(:,i) = (A*V(:,i)) -(alph(i)*V(:,i)) - (bet(i)*V(:,i-1)); 
    
    end
%Calculating the tridiagonal matrix

% Setting the initial values
T(1,1) = alph(1);
T(1,2) = bet(2);
T(m,m) = alph(m);
T(m,m-1) = bet(m);
T(m-1,m) = bet(m);

%Setting the non terminal values
for i=2:m-1
    T(i,i) = alph(i);
    T(i,i-1) = bet(i);
    T(i,i+1) = bet(i+1);
    
end
    [E1,~] = eigs(T,1,'sr'); % smallest real eigenvector
    E1 = V*E1;
    V=E1;
    E = diag((V')*A*V);
else
    [V,E] = eigs(A,1,'sr');
    T=A;
end
    

    
end
