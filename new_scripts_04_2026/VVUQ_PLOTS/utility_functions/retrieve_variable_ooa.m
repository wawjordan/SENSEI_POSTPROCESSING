function variable = retrieve_variable_ooa( ALL_DATA, folder_num, geom, prim, var )
    error = abs(retrieve_variable_error( ALL_DATA, folder_num, geom, prim, var ));
    N = retrieve_variable( ALL_DATA, folder_num, geom, '', 'N' );
    variable = nan*error;
    r_fac = 2;
    for i = 2:numel(N)
        variable(:,i) = log(error(:,i-1) ./ error(:,i))./log(r_fac);
    end
end