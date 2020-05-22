function [A] = canon(M,j) %converts MPS in site canonical form
%input check
if ((size(M,2) < j))
    error('ERROR: the site specified is beyond the dimension of MPS.');
end
N = size(M,2);
% end points
if ( j ==1)
    
    for i = (N:-1:(j+1))
        
    
    
        T = M{i};
        T = reshape(T, [size(T,1), size(T,2)*size(T,3)]);
        [U,S,V] = svd(T,'econ');
        M{i} = reshape(V',[size(V,2), size(M{i},2), size(V,1)/size(M{i},2)]);

        if (i== j+1)
            
            M{i-1} = ncon({(M{i-1}),U*S},{[-1 -2 1],[1 -3]});
        end
    end
else

     %convert a tensor in left canonical form

for i = (1:j-1)
    T = M{i}; % tensor at ith site 
    T = reshape(T, [size(T,1)*size(T,2), size(T,3)]); % Convert the rank 3 tensor into a matrix to perform SVD
    [U,S,V] = svd(T,'econ'); % econ helps in memory 
    M{i} = reshape(U,[size(U,1)/size(M{i},2), size(M{i},2), size(U,2)]); % Replacing the tensor with the U 
    if (i== j-1)
         
        M{i+1} = ncon({(S*V'),M{i+1}},{[-1 1],[1 -2 -3]});
        
    end
end







for i = (N:-1:(j+1))
    
    
    T = M{i};
    T = reshape(T, [size(T,1), size(T,2)*size(T,3)]);
    [U,S,V] = svd(T,'econ');
    M{i} = reshape(V',[size(V,2), size(M{i},2), size(V,1)/size(M{i},2)]);

if (i== j+1)
    M{i-1} = ncon({(M{i-1}),U*S},{[-1 -2 1],[1 -3]});
end
end



end
A=M;
end