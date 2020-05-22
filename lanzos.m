%Implement lanczos algorithm for a matrix 
%Inputs a matrix(n*n) and number of iterations and the initial seed
%Outputs V(n*m) with orthonormal columns and betas and alphas
m =5;
A=diag([2 1 -1 -2 5]);
T=zeros(m,m);
r = (rand(5,1));
alph = zeros(1,m);
bet = zeros(1,m);
[s1,~] = size(A);
V = zeros(s1,m);
r= r/norm(r);
W = zeros(s1,m);    
alph(1) = (r.')*(A.')*(r);
V(:,1) = r;
W(:,1) = A*r - alph(1)*r;
for i= 2:m
    bet(i) = norm(W(:,i-1));
    V(:,i) = (W(:,i-1))/ bet(i); 
    alph(i) = (V(:,i).')*(A.')*(V(:,i));
    W(:,i) = (A*V(:,i)) -(alph(i)*V(:,i)) - (bet(i)*V(:,i-1)); 
    
end   

T(1,1) = alph(1);
T(1,2) = bet(2);
T(m,m) = alph(m);
T(m,m-1) = bet(m);
for i=2:m-1
    T(i,i) = alph(i);
    T(i,i-1) = bet(i);
    T(i,i+1) = bet(i+1);
end
    
e = eig(T);







