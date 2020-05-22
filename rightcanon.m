function [B] = rightcanon(M)

N=size(M,2);

for i = (N:-1:1)
    
    
    T = M{i};
    T = reshape(T, [size(T,1), size(T,2)*size(T,3)]);
    [U,S,V] = svd(T,'econ');
    M{i} = reshape(V',[size(V,2), size(M{i},2), size(V,1)/size(M{i},2)]);

if i > 1
    M{i-1} = ncon({(M{i-1}),U*S},{[-1 -2 1],[1 -3]});
end
if i == 1
    M{i} = reshape(S*V',[size(V,2), size(M{i},2), size(V,1)/size(M{i},2)]);
end
end


B=M;
end