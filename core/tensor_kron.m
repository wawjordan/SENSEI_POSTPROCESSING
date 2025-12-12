function out = tensor_kron(A,B)
n_dim = numel(size(A));

order = fliplr( reshape(2*n_dim:-1:1,n_dim,2).' );

out = reshape( permute(tensorprod(A,B),order(:)), size(A).*size(B) );
end