function [A] = leftcanon(M) %convert a tensor in left canonical form
N = size(M,2);
for i = (1:N)
    T = M{i}; % tensor at ith site 
    T = reshape(T, [size(T,1)*size(T,2), size(T,3)]); % Convert the rank 3 tensor into a matrix to perform SVD
    [U,S,V] = svd(T,'econ'); % econ helps in memory 
    M{i} = reshape(U,[size(U,1)/size(M{i},2), size(M{i},2), size(U,2)]); % Replacing the tensor with the U 
    if i<N
         
        M{i+1} = ncon({(S*V'),M{i+1}},{[-1 1],[1 -2 -3]});
    end
    if i==N
         
        M{i} = reshape(U*S,[size(U,1)/size(M{i},2), size(M{i},2), size(U,2)]);
    end
end
A = M;
end