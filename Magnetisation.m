function [O] = Magnetisation(A,a)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
t = size(A,2);
if(a<1 || a>t)
    error('ERR: Outof bound values');
end

s = (size(A{1},2)-1)/2;
[Sz,~,~,~] = Spin(s);

if(a~=1)
   A1= ncon({A{1},conj(A{1})},{[1,2,-1],[1,2,-2]});
else
   A1= ncon({A{1},Sz,conj(A{1})},{[1,2,-1],[2,3],[1,3,-2]}); 
end


        
     

for i=(2:a-1)
    A1= ncon({A1,A{i},conj(A{i})},{[1,2],[1,3,-1],[2,3,-2]});
end

for i=(a:t)
    if (i==a && a~=1)
        A1= ncon({A1,A{i},conj(A{i}),Sz},{[1,2],[1,3,-1],[2,4,-2],[3,4]});
    elseif(i~=a)
        A1= ncon({A1,A{i},conj(A{i})},{[1,2],[1,3,-1],[2,3,-2]});        
    end
end

O=A1;

end



