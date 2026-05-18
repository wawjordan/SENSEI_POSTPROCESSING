function variable = retrieve_variable( ALL_DATA, folder_num, geom, prim, var )
if strcmp(var,'N')
    variable = [ALL_DATA.DATA(folder_num).(geom).F(:).N];
    variable = variable/ALL_DATA.DATA(folder_num).n_skip;
elseif strcmp(var,'XC')
    variable = {ALL_DATA.DATA(folder_num).(geom).F(:).XC};
elseif strcmp(var,'xc')||strcmp(var,'yc')||...
       strcmp(var,'x') ||strcmp(var,'y')||...
       strcmp(var,'x0')||strcmp(var,'y0')
    N = numel(ALL_DATA.DATA(folder_num).(geom).G);
    variable = cell(N,1);
    for i = 1:N
        variable{i} = ALL_DATA.DATA(folder_num).(geom).G(i).(var);
    end
elseif strcmpi(var,'theta_interval')
    airfoil = make_airfoil(ALL_DATA);
    N = numel(ALL_DATA.DATA(folder_num).(geom).G);
    variable = cell(N,1);
    for i = 1:N
        X = ALL_DATA.DATA(folder_num).(geom).G(i).x;
        Y = ALL_DATA.DATA(folder_num).(geom).G(i).y;
        variable{i} = airfoil.get_theta_from_coords_piecewise_constant(X,Y,true);
    end
elseif strcmpi(var,'theta_center')
    airfoil = make_airfoil(ALL_DATA);
    N = numel(ALL_DATA.DATA(folder_num).(geom).G);
    variable = cell(N,1);
    for i = 1:N
        X = ALL_DATA.DATA(folder_num).(geom).G(i).xc;
        Y = ALL_DATA.DATA(folder_num).(geom).G(i).yc;
        tmp = airfoil.get_theta_from_coords(X,Y,true);
        variable{i} = tmp(1:2:end);
    end
elseif strcmp(var,'CP')||strcmp(var,'P')
    N = numel(ALL_DATA.DATA(folder_num).(geom).F);
    variable = cell(N,1);
    for i = 1:N
        variable{i} = ALL_DATA.DATA(folder_num).(geom).F(i).([prim,'_',var]);
    end
else
    variable = [ALL_DATA.DATA(folder_num).(geom).H(:).([prim,'_',var])];
end
end