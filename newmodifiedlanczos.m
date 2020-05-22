% m is the number of lanzos vectors considered for the diagonalisation
% if ((size(A,1) ~= size(A,2)))
%     error('ERROR: A is not square matrix.');
% elseif size(A,1) ~= numel(r)
%     error('ERROR: r is not compatible with A.');
% end
function [y,E] = newmodifiedlanczos(r,L,R,MO)
global Ld
l= size(L,1);
d = size(MO,2);
r2 = size(R,1);
if( l~=0 && r2~=0)
    s1 = l*r2*d;
elseif(l==0)
    s1 =d*r2;
elseif(r2==0)
    s1= l*d;  
end
%A= (A +A')/2; % To reduce numerical noise
m= Ld;
if ((s1<=m))
    m=s1;
end

    
    
r= r/norm(r);% normalising the seed vector
%[s1,~] = size(A);

Ar = product(r,MO,L,R); %try to make the lanczos efficient

V = zeros(s1,m);
alph1 = (Ar')*r;
%alph1 = (r')*(A')*r; % Setting the first value
bet2 =norm((Ar) - (alph1*r));
%bet2 =norm((A*r) - (alph1*r));
alph2 = product((Ar-alph1*r),MO,L,R);
alph2 =((alph2')*((Ar - alph1*r))/(bet2*bet2));
%alph2 =(((A*r - alph1*r)')*(A')*((A*r - alph1*r))/(bet2*bet2));


% initialising the tri matrix
T = zeros(m,m);
T(1,1) = alph1;
T(1,2) = bet2;
T(2,1) = bet2;
T(2,2) = alph2;

%Initial vectors
v1 =r;
v2 =  (Ar - alph1*r)/bet2;
%v2 =  (A*r - alph1*r)/bet2;
V(:,1) = v1;
V(:,2) = v2;



for i= 1:m
    if( (i==1) || (i==2))
        vpp = v1; % vector from previous to previous iteration
        vp = v2; % vector from previous iteration
    elseif i ==3
        Avp = product(vp,MO,L,R);
        T(i-1,i) = norm( (Avp) - ((T(i-1,i-1))*vp) - ((T(i-2,i-1))*vpp) ); % setting the value of beta
        % T(i-1,i) = norm( (A*vp) - ((T(i-1,i-1))*vp) - ((T(i-2,i-1))*vpp) ); % setting the value of beta
      
        T(i,i-1) = norm( (Avp) - ((T(i-1,i-1))*vp) - ((T(i-2,i-1))*vpp) );% setting the value of beta
        %T(i,i-1) = norm( (A*vp) - ((T(i-1,i-1))*vp) - ((T(i-2,i-1))*vpp) );% setting the value of beta
        %bet(i) = norm((A*vp) - ((T(i-1,i-1))*vp) - ((T(i-1,i))*vpp));
        
        vc = (Avp - alph2*vp - bet2*v1)/(T(i-1,i));
        %vc = (A*vp - alph2*vp - bet2*v1)/(T(i-1,i));
        
        V(:,i) = vc;
        %alph(i) = (vc')*(A')*(vc);
        
        T(i,i) = (product(vc,MO,L,R)')*(vc);
        %T(i,i) = (vc')*(A')*(vc);
    else 
        vpp = vp; % modify the previous vectors according to the iteration
        vp = vc; % modify the previous vectors according to the iteration
        Avp = product(vp,MO,L,R);
        
        T(i-1,i) = norm((Avp) - ((T(i-1,i-1))*vp) - ((T(i-2,i-1))*vpp));
        %T(i-1,i) = norm((A*vp) - ((T(i-1,i-1))*vp) - ((T(i-2,i-1))*vpp));
        
        T(i,i-1) = norm((Avp) - ((T(i-1,i-1))*vp) - ((T(i-2,i-1))*vpp));
        %T(i,i-1) = norm((A*vp) - ((T(i-1,i-1))*vp) - ((T(i-2,i-1))*vpp));
        
        vc = (Avp - (T(i-1,i-1))*vp - (T(i-2,i-1))*vpp)/(T(i-1,i));
        %vc = (A*vp - (T(i-1,i-1))*vp - (T(i-2,i-1))*vpp)/(T(i-1,i));
        V(:,i) = vc;
        %alph(i) = (vc')*(A')*(vc);
        
        T(i,i) = (product(vc,MO,L,R)')*(vc); % setting the value of alphas
        %T(i,i) = (vc')*(A')*(vc); % setting the value of alphas
    
    end
   
    
end
     T = 0.5*(T+T'); % reduce numerical noise
    [y,E] = eigs(T,1,'sr'); % 'sr': smallest real
    y = V*y; 
    
    %Ay = product(y,MO,L,R);
    
    %E = diag((y')*Ay); % better approximation of energy
    %E = diag((y')*A*y); % better approximation of energy
end
