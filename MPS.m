N = 10; % number of sites
d = 2; % local space dimension
 global D % bond dimension when truncation happens (bond dimension of the MPS)
 M = cell(1,N); % representation of a MPS; M{n} is the tensor at site n
% s = floor(log(D)/log(d)); %site at r =rand(1,2,which truncation happens
%  for j=(1:N/2)
%      if j == 1
%          M{j} = rand(1,d,d); % left end; left leg is of size 1
%          M{N} = rand(d,d,1);
%      elseif (d^j<=D) 
%          M{j} = rand((d^(j-1)),d,(d^j)); % bond dimension before truncation of left tensors
%          M{N-j+1} = rand((d^j),d,(d^(j-1)));
%      elseif(d^j>D && (d^(j-1)<=D)) 
%          M{j} = rand((d^(j-1)),d,D); % bond dimension after truncation
%          M{N-j+1} = rand(D,d,(d^(j-1)));
%      else
%          M{j} = rand(D,d,D); % bond dimension after truncation
%          M{N-j+1} = rand(D,d,D);
%      end
%  end

        
        
    

    
    

 for i = (1:N)
  %assign individual tensors
 % leg order: left, bottom, right
 if i == 1
% % left end; left leg is of size 1
 M{i} = rand(1,d,D);
 elseif i == N
% % right end; right leg is of size 1
 M{i} = rand(D,d,1);
 else
 M{i} = rand(D,d,D);
 end
 end