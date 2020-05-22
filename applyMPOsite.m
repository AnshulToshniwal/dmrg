function [Aeff,Eeff] = applyMPOsite(Hleft,Hloc,Hright)

if isempty(Hleft) && isempty(Hright)
    error('Error: Hleft and Hright cannot be empty at the same time.');
elseif ~isempty(Hleft) && (size(Hleft,1) ~= size(Hleft,3))
    error('Error: the first and third legs of Hleft should be of the same dimension.');
elseif ~isempty(Hright) && (size(Hright,1) ~= size(Hright,3))
    error('Error: the first and third legs of Hright should be of the same dimension.');
end
% % %

% if ~isempty(Hleft)
%     Heff = ncon({Hleft,Hloc},{[-1 1 -2],[1 -3 -4 -5]}); 
% else
%     Heff = reshape(Hloc,[1 size(Hloc)]); % added a  dummy leg, for convenience
% end
% if ~isempty(Hright)
%     Heff = ncon({Heff,Hright},{[-1 -2 -3 1 -4],[-5 1 -6]}); 
% else
%     Heff = permute(Heff,[1 2 3 5 4]); % permute the 4th leg (right leg of Hloc) to the last
% end
% 
% t = zeros(1,3);
% for i = (1:numel(t))
%     t(i) = size(Heff,2*i);
% end
% Heff = permute(Heff,[(1:2:6),(2:2:6)]);
% Heff = reshape(Heff,[1 1]*prod(t));
% s1 = size(Heff,1);

l= size(Hleft,1);
d = size(Hloc,2);
r2 = size(Hright,1);
if( l~=0 && r2~=0)
    s1 = l*r2*d;
    r = rand(s1,1);
   [Aeff,Eeff] = newmodifiedlanczos(r,Hleft,Hright,Hloc); % hope it doesnt make a mess
   Aeff = reshape(Aeff,l,d,r2);
elseif(l==0)
    s1 =d*r2;
    r = rand(s1,1);
    [Aeff,Eeff] = newmodifiedlanczos(r,Hleft,Hright,Hloc); % hope it doesnt make a mess
    Aeff = reshape(Aeff,1,d,r2);
elseif(r2==0)
    s1= l*d; 
    r = rand(s1,1);
    [Aeff,Eeff] = newmodifiedlanczos(r,Hleft,Hright,Hloc); % hope it doesnt make a mess
    Aeff = reshape(Aeff,l,d,1);
end





 %Replace by Lanczos 
%r = rand(s1,1);
%[Aeff,Eeff] = modifiedlanczos(Heff,r); % hope it doesnt make a mess
%[Aeff,Eeff] = eigs(Heff,1,'sr');

%Aeff = reshape(Aeff,t);


%need to reshape the matrix accordin to local space dimension
%alpha = size(Hleft,1);
%beta = size(Hright,1);

%s =(size(Hloc,2));
%V = reshape(V1,[alpha*(s),beta*(s)]);

end