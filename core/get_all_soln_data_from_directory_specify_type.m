function [S,file_names] = get_all_soln_data_from_directory_specify_type(folder,pat1,pat2,read_again)
% pat1 = '*-soln.mat';
% pat2 = '*-soln.dat';
S = struct();
TMP1 = dir(fullfile(folder,'**\',pat1));
N_files1 = length(TMP1);
folders1 = cell(N_files1,1);

file_names = cell(N_files1,1);

if (nargin==3)
    read_again=false;
end

if ~read_again
    if ~isempty(TMP1)
        % get the corresponding folders and soln (.mat) data
        for i = 1:N_files1
            folders1{i} = replace(TMP1(i).folder,{folder,pat1},{'',''});
            S(i).N  = cellfun(@str2double,regexp(folders1{i},'\d*','match'));
            S(i).N  = sqrt(prod(S(i).N));
            file_name = fullfile(TMP1(i).folder,TMP1(i).name);
            load(file_name,'DATA');
            S(i).DATA = DATA;
        end
    end
end


TMP2 = dir(fullfile(folder,'**\',pat2));
N_files2 = length(TMP2);
folders2 = cell(N_files2,1);
% get the corresponding folders and soln (.dat) data (if not already loaded)
for i = 1:N_files2
    folders2{i} = replace(TMP2(i).folder,{folder,pat2},{'',''});
    
    file_name = fullfile(TMP2(i).folder,TMP2(i).name);
    file_name2 = replace(file_name, replace(pat2,'*',''), replace(pat1,'*','') );

    file_names{i} = file_name;
    % don't read in the same data again
    if any(cellfun(@(s2)strcmp(folders2{i},s2),folders1))
        continue
    end
    S(i).N  = cellfun(@str2double,regexp(folders2{i},'\d*','match'));
    S(i).N  = sqrt(prod(S(i).N));
    % DATA = tec2mat_structured_condensed(file_name);
    DATA = get_SENSEI_tecplot_soln_data(file_name);
    save(file_name2,'DATA');
    S(i).DATA = DATA;
end

end