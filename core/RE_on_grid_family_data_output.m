function [ S, S_ERR_exact, S_RE_out, S_ERR_out, S_RE_ERR_out ] = RE_on_grid_family_data_output(foldername,p,var_list,read_again,has_exact)
pat1 = '*-soln.dat';
pat2 = '*-soln_RE.dat';
pat3 = '*-soln_RE_ERR.dat';
pat4 = '*-soln_RE_ERR_EXACT.dat';
% N_files1 = length(TMP1);
% folders1 = cell(N_files1,1);


fprintf('reading data from directory\n');
[S,file_names] = get_all_soln_data_from_directory(foldername,read_again);

file_names = file_names(2:end);

fprintf('performing Richardson Extrapolation\n');
[ S_RE_out, S_ERR_out ] = RE_on_grid_family_data(S,p,var_list);
S_RE_ERR_out = S_ERR_out;

if (has_exact)
    [S_ERR_exact,~] = get_all_soln_data_from_directory_specify_type(foldername,'*-DE-primitive_TEST.mat','*-DE-primitive_TEST.dat',read_again);
    for i = 2:numel(S_ERR_exact)
        for v = 1:numel(var_list)
            if (isfield(S_ERR_exact(i).DATA(1).ZONE,var_list{v}))
                S_RE_ERR_out(i-1).DATA.ZONE.(var_list{v}) = S_ERR_out(i-1).DATA.ZONE.(var_list{v}) - S_ERR_exact(i).DATA(1).ZONE.(var_list{v});
            end
        end
    end
end

fprintf('writing out data\n');
for i = 1:numel(file_names)
    filename = replace(file_names{i}, replace(pat1,'*',''), replace(pat2,'*','') );
    write_tecplot_zones_from_DATA(filename,S_RE_out(i).DATA)

    filename = replace(file_names{i}, replace(pat1,'*',''), replace(pat3,'*','') );
    write_tecplot_zones_from_DATA(filename,S_ERR_out(i).DATA)

    if ( has_exact)
        filename = replace(file_names{i}, replace(pat1,'*',''), replace(pat4,'*','') );
        write_tecplot_zones_from_DATA(filename,S_RE_ERR_out(i).DATA)
    end
end

end