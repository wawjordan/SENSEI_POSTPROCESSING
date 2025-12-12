function [ S_RE_out, S_ERR_out ] = RE_on_grid_family_data(S,p,var_list)

Nzones = numel(S);
n_vars = numel(var_list);

%                                             1   -->  Nzones
% assumes finer grids in ascending order ( coarse -->   fine    )


% preallocate
S_RE_out  = S(2:Nzones);
S_ERR_out = S(2:Nzones);

for n = 2:Nzones
    
% [RE_vals,ERR_vals] = perform_Richardson_Extrapolation_known_p( fine_zone, coarse_zone, p, var_list );
[RE_vals,ERR_vals] = perform_Richardson_Extrapolation_known_p( S(n).DATA.ZONE, S(n-1).DATA.ZONE, p, var_list );

tmp_zone1 = coarse_data_to_fine_grid_inject(S(n).DATA.ZONE,RE_vals,var_list);
tmp_zone2 = coarse_data_to_fine_grid_inject(S(n).DATA.ZONE,ERR_vals,var_list);
for v = 1:n_vars
   S_RE_out(n-1).DATA.ZONE.(var_list{v}) = tmp_zone1.(var_list{v});
   S_ERR_out(n-1).DATA.ZONE.(var_list{v}) = tmp_zone2.(var_list{v});
end

end