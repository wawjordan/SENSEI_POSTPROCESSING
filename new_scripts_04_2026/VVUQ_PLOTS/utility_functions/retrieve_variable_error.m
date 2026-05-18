function variable = retrieve_variable_error( ALL_DATA, folder_num, geom, prim, var )
    variable = retrieve_variable( ALL_DATA, folder_num, geom, prim, var );
    airfoil = make_airfoil(ALL_DATA);
    variable = variable-airfoil.(var);
end