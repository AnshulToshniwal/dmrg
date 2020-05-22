function [A] = product(r3,MO,Ls,Rs)% function to evaluate the product H*v
 MP = reshape(r3,[size(Ls,3),size(MO,2),size(Rs,3)]);
if (~isempty(Ls) && ~isempty(Rs)) 
 

  
  %A= ncon({Ls,MP,MO,Rs},{[1,2,-1],[1,3,4],[2,3,5,-2],[-3,5,4]},[2 4 5 1 3]);% came pretty close!
  A= ncon({Ls,MP,MO,Rs},{[-1,1,2],[2,3,4],[1,-2,5,3],[-3,5,4]},[2 4 3 1 5]); % yeah BITCH!!!
  A = reshape(A,(size(A,1)*size(A,2)*size(A,3)),1);
 
elseif(isempty(Rs))
    A= ncon({Ls,MP,MO},{[-1,1,2],[2,3,-2],[1,-3,-4,3]},[2 3 1]);
    A = reshape(A,(size(A,1)*size(A,2)*size(A,3)),1);
else
    A= ncon({MP,MO,Rs},{[-1,1,3],[-4,-2,2,1],[-3,2,3]},[2 3 1]);
    A = reshape(A,(size(A,1)*size(A,2)*size(A,3)),1);
 
end
end 