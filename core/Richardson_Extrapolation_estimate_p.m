function [p_avg,p_vals_zone] = Richardson_Extrapolation_estimate_p( fine_zone, med_zone, coarse_zone, var_list, p_formal )

if nargin<4
    p_formal = 2;
end

min_p = 0.5;

n_vars = numel(var_list);
p_avg = zeros(n_vars,1);

[fine_on_coarse,~] = fine_data_to_coarse_grid(fine_zone,coarse_zone,var_list);
[med_on_coarse,rh] = fine_data_to_coarse_grid(med_zone,coarse_zone,var_list);

r = max(rh); % ?

p_vals_zone = med_on_coarse;

% phat = real(log((Grid3Data(:,:,vv) - Grid2Data(:,:,vv)) ./ ...
%                     (Grid2Data(:,:,vv) - Grid1Data(:,:,vv))) / log(r));
%     phat = min(max(phat, 0.5), 2);
%     pglb(vv) = mean(phat(:),'omitnan');

omit = @(a) (~isreal(a))&(a<=0);
vol = coarse_zone.('VOL');
for v = 1:n_vars
    p_tmp  = log(  ( coarse_zone.(var_list{v}) - med_on_coarse.(var_list{v}) ) ...
                 ./( med_on_coarse.(var_list{v}) - fine_on_coarse.(var_list{v}) ) ) / log( r );
    mask = omit(p_tmp);
    p_tmp(mask) = nan;
    p_tmp(~mask) = min( max( p_tmp(~mask), min_p), p_formal);
    mask2 = isnan(p_tmp);
    p_avg(v) = sum( p_tmp(~mask2).*vol(~mask2),"all")/ sum(vol(~mask2),"all");
    p_tmp(mask2) = 0.0;
    p_vals_zone.(var_list{v}) = p_tmp;
end

end