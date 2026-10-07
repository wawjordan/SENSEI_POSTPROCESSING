function [E,nx,ni] = concatenate_error_norms_from_struct(DE_struct,var,norm)
n_grids = numel(DE_struct);
n_iters = [DE_struct.N_iterations];
ni = 0:max(n_iters);
nx = [DE_struct.N2];
E = nan(n_grids,max(n_iters)+1);
for i = 1:n_grids
    E(i,1:n_iters(i)+1) = squeeze(DE_struct(i).E(var,norm,:));
end
end

