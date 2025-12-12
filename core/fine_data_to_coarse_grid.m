function [coarse_zone_out,refine_ratio] = fine_data_to_coarse_grid(fine_zone,coarse_zone,var_list)
n_vars = numel(var_list);
inputs = cell(1,2*n_vars);

cnt = 0;
for v = 1:n_vars
    cnt = cnt + 1;
    % clean up names
    inputs{cnt} = var_list{v};
    cnt = cnt + 1;
    inputs{cnt} = coarse_zone.(var_list{v});
end

coarse_zone_out = struct(inputs{:});


refine_ratio = size(fine_zone.('VOL'))./size(coarse_zone.('VOL'));

for v = 1:numel(var_list)
    coarse_zone_out.(var_list{v}) = get_coarse_cell_avg( fine_zone.('VOL'),fine_zone.(var_list{v}),refine_ratio);
end

end


function agglom = get_coarse_cell_avg(fine_volume,fine_data,n_skip)
sz = size(fine_volume);
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
ssz = sz(:)./n_skip_(:);

dim_dist = cell(n_dim,1);
for i = 1:n_dim
    dim_dist{i} = n_skip_(i)*ones(ssz(i),1);
end

agglom = cellfun( @(A) sum(A,"all"), mat2cell( fine_data.*fine_volume, dim_dist{:} ) ) ...
       ./cellfun( @(A) sum(A,"all"), mat2cell( fine_volume, dim_dist{:} ) );

end