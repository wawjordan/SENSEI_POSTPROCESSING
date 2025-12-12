function [fine_zone_out,refine_ratio] = coarse_data_to_fine_grid_inject(fine_zone,coarse_zone,var_list)
n_vars = numel(var_list);
inputs = cell(1,2*n_vars);

cnt = 0;
for v = 1:n_vars
    cnt = cnt + 1;
    % clean up names
    inputs{cnt} = var_list{v};
    cnt = cnt + 1;
    inputs{cnt} = fine_zone.(var_list{v});
end

fine_zone_out = struct(inputs{:});

refine_ratio = size(fine_zone.(var_list{1}))./size(coarse_zone.(var_list{1}));

for v = 1:numel(var_list)
    fine_zone_out.(var_list{v}) = get_coarse_cell_avg( coarse_zone.(var_list{v}),refine_ratio);
end

end


function finer = get_coarse_cell_avg(coarse_data,n_skip)
sz = size(coarse_data);
n_dim = numel(sz);
if ~isscalar(n_skip)
    n_skip_ = ones(n_dim,1);
    nnskip = min(n_dim,numel(n_skip));
    n_skip_(1:nnskip) = n_skip(1:nnskip);
else
    n_skip_ = n_skip*ones(n_dim,1);
end

if any(mod(sz,n_skip_)~=0)
    error('size of array must be divisible by n_skip')
end

dim = num2cell(n_skip_);
% really need to check that this is correct ...
finer = tensor_kron( coarse_data, ones(dim{:}) );

end

function out = tensor_kron(A,B)
n_dim = numel(size(A));

order = fliplr( reshape(2*n_dim:-1:1,n_dim,2).' );

out = reshape( permute(tensorprod(A,B),order(:)), size(A).*size(B) );
end