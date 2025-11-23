function S = get_soln_data_from_directory(folder)
pat1 = '*-soln.mat';
pat2 = '*-soln.dat';
S = struct();
TMP1 = dir(fullfile(folder,'**\',pat1));
N_files1 = length(TMP1);
folders1 = cell(N_files1,1);

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

TMP2 = dir(fullfile(folder,'**\',pat2));
N_files2 = length(TMP2);
folders2 = cell(N_files2,1);
% get the corresponding folders and soln (.dat) data (if not already loaded)
for i = 1:N_files2
    folders2{i} = replace(TMP2(i).folder,{folder,pat2},{'',''});
    % don't read in the same data again
    if any(cellfun(@(s2)strcmp(folders2{i},s2),folders1))
        continue
    end
    S(i).N  = cellfun(@str2double,regexp(folders2{i},'\d*','match'));
    S(i).N  = sqrt(prod(S(i).N));
    file_name = fullfile(TMP2(i).folder,TMP2(i).name);
    DATA = tec2mat_structured_condensed(file_name);
    file_name2 = replace(file_name, replace(pat2,'*',''), replace(pat1,'*','') );
    save(file_name2,'DATA');
    S(i).DATA = DATA;
end

end