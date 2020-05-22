
%Outputs the tridiagonal matrix
function [y,E] = modifiedlanczos(A,r)
%m is the number of lanzos vectors considered for the diagonalisation
if ((size(A,1) ~= size(A,2)))
    error('ERROR: A is not square matrix.');
elseif size(A,1) ~= numel(r)
    error('ERROR: r is not compatible with A.');
end
global Ld
[s1,~] = size(A);
A= (A +A')/2; % To reduce numerical noise
m= Ld;
if ((s1<=m))
    [y,E] = eigs(A,1,'sr');
    
    
else

    
    
r= r/norm(r);% normalising the seed vector




V = zeros(s1,m);
alph1 = (r')*(A')*r; % Setting the first value
bet2 =norm((A*r) - (alph1*r));
alph2 =(((A*r - alph1*r)')*(A')*((A*r - alph1*r))/(bet2*bet2));


% initialising the tri matrix
T = zeros(m,m);
T(1,1) = alph1;
T(1,2) = bet2;
T(2,1) = bet2;
T(2,2) = alph2;

%Initial vectors
v1 =r;
v2 =  (A*r - alph1*r)/bet2;
V(:,1) = v1;
V(:,2) = v2;



for i= 1:m
    if( (i==1) || (i==2))
        vpp = v1; % vector from previous to previous iteration
        vp = v2; % vector from previous iteration
    elseif i ==3
        
        
        T(i-1,i) = norm( (A*vp) - ((T(i-1,i-1))*vp) - ((T(i-2,i-1))*vpp) ); % setting the value of beta
        T(i,i-1) = norm( (A*vp) - ((T(i-1,i-1))*vp) - ((T(i-2,i-1))*vpp) );% setting the value of beta
        %bet(i) = norm((A*vp) - ((T(i-1,i-1))*vp) - ((T(i-1,i))*vpp));
        
        
        vc = (A*vp - alph2*vp - bet2*v1)/(T(i-1,i));
        
        V(:,i) = vc;
        %alph(i) = (vc')*(A')*(vc);
        
        
        T(i,i) = (vc')*(A')*(vc);
    else 
        vpp = vp; % modify the previous vectors according to the iteration
        vp = vc; % modify the previous vectors according to the iteration
        
        
        
        T(i-1,i) = norm((A*vp) - ((T(i-1,i-1))*vp) - ((T(i-2,i-1))*vpp));
        T(i,i-1) = norm((A*vp) - ((T(i-1,i-1))*vp) - ((T(i-2,i-1))*vpp));
        
        
        vc = (A*vp - (T(i-1,i-1))*vp - (T(i-2,i-1))*vpp)/(T(i-1,i));
        V(:,i) = vc;
        %alph(i) = (vc')*(A')*(vc);
        
        
        T(i,i) = (vc')*(A')*(vc); % setting the value of alphas
    
    end
   
    
end
     T = 0.5*(T+T'); % reduce numerical noise
     [y,E] = eigs(T,1,'sr'); % 'sr': smallest real
     y = V*y; 
    % E = diag((y')*A*y); % better approximation of energy
end
end