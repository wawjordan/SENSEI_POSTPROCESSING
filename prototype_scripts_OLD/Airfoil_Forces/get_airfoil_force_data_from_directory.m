function DATA = get_airfoil_force_data_from_directory(folder,alpha)
% iterative correction (no offsets)
history_pat = { '*-001-exact-force_history.dat',  ...% from exact cell avgs
                '*-001-primal-force_history.dat', ...% from primal cell avgs
                '*-001-ete-force_history.dat',    ...% from ete corrected cell avgs
                '*-001-ic-force_history.dat'};       % from iterative corrected cell avgs
surface_pat = { '*-exact-surface_forces.dat',     ...%
                '*-primal-surface_forces.dat',    ...%
                '*-ete-surface_forces.dat',       ...%
                '*-ic-surface_forces.dat'};
DATA = struct();
for i = 1:4
    DATA(i).F = read_surface_forces_from_folder(folder,surface_pat{i});
    DATA(i).H = read_force_histories_from_folder(folder,history_pat{i},alpha);
end
end

function F = read_surface_forces_from_folder(folder,pattern)
F = struct();
TMP = dir(fullfile(folder,'**\',pattern));
N_files = length(TMP);
if (N_files==0)
    return
end
for i = 1:N_files
    file_name = fullfile(TMP(i).folder,TMP(i).name);
    [F(i).XC,F(i).P] = get_airfoil_surface_pressure_data_from_file(file_name);
    str = replace(TMP(i).folder,folder,'');
    F(i).N = cellfun(@str2double,regexp(str,'\d*','match'));
    F(i).file_name = file_name;
    F(i).name = pattern;
end
end

function H = read_force_histories_from_folder(folder,pattern,alpha)
H = struct();
TMP = dir(fullfile(folder,'**\',pattern));
N_files = length(TMP);
if (N_files==0)
    return
end
for i = 1:N_files
    file_name = fullfile(TMP(i).folder,TMP(i).name);
    [H(i).CX,H(i).CY,H(i).CL,H(i).CD] = get_force_history_data_from_file(file_name,alpha);
    str = replace(TMP(i).folder,folder,'');
    H(i).N = cellfun(@str2double,regexp(str,'\d*','match'));
    H(i).file_name = file_name;
    H(i).name = pattern;
end
end