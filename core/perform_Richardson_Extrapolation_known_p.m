function [RE_vals,ERR_vals] = perform_Richardson_Extrapolation_known_p( fine_zone, coarse_zone, p, var_list )
n_vars = numel(var_list);
p_hat = zeros(n_vars,1);
if isscalar(p)
    p_hat = 0*p_hat + p;
else
    if (numel(p)==numel(p_hat))
        p_hat = p;
    else
        error('must provide a scalar p, or a list with the samne number of elements as var_list.')
    end
end

[fine_data_on_coarse_grid,r] = fine_data_to_coarse_grid(fine_zone,coarse_zone,var_list);

r = max(r); % ?

RE_vals = fine_data_on_coarse_grid;
ERR_vals = fine_data_on_coarse_grid;

% RichExtrapSol(:,:,vv) = Grid1Data(:,:,vv) + ...
%                             (Grid1Data(:,:,vv) - Grid2Data(:,:,vv)) / (r^p - 1);
% 
%     DiscError(:,:,vv) = RichExtrapSol(:,:,vv) - Grid1Data(:,:,vv);
for v = 1:n_vars
    ERR_vals.(var_list{v}) = -( fine_data_on_coarse_grid.(var_list{v}) ...
                             - coarse_zone.(var_list{v}) ) / (r^p_hat(v)-1);
    RE_vals.(var_list{v}) =  fine_data_on_coarse_grid.(var_list{v}) ...
                          + ERR_vals.(var_list{v});
end

end