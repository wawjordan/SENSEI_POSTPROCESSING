function [N,E,OOA] = gather_error_plot_data(ALL_DATA,folders,geom,prim,var,ics)
N = {};
E = {};
OOA = {};
for i = 1:numel(folders)
    tmp = retrieve_variable( ALL_DATA, folders(i), geom{i}, '', 'N' );
    N = [N,tmp];
    tmp = abs(retrieve_variable_error( ALL_DATA, folders(i), geom{i}, prim{i}, var{i} ));
    tmp = tmp(ics{i}(:)+1,:);
    E = [E,tmp];
    tmp = retrieve_variable_ooa( ALL_DATA, folders(i), geom{i}, prim{i}, var{i} );
    tmp = tmp(ics{i}(:)+1,:);
    OOA = [OOA,tmp];
end
end