function [ S, S_RE_out, S_ERR_out ] = RE_on_grid_family_data_output(foldername,p,var_list,read_again)
pat1 = '*-soln.dat';
pat2 = '*-soln_RE.dat';
pat3 = '*-soln_RE_ERR.dat';
% N_files1 = length(TMP1);
% folders1 = cell(N_files1,1);


fprintf('reading data from directory\n');
[S,file_names] = get_all_soln_data_from_directory(foldername,read_again);
file_names = file_names(2:end);

fprintf('performing Richardson Extrapolation\n');
[ S_RE_out, S_ERR_out ] = RE_on_grid_family_data(S,p,var_list);

fprintf('writing out data\n');
for i = 1:numel(file_names)
    filename = replace(file_names{i}, replace(pat1,'*',''), replace(pat2,'*','') );
    write_tecplot_zones_from_DATA(filename,S_RE_out(i).DATA)

    filename = replace(file_names{i}, replace(pat1,'*',''), replace(pat3,'*','') );
    write_tecplot_zones_from_DATA(filename,S_ERR_out(i).DATA)
end

end
% function S = get_all_soln_data_from_directory(folder,read_again)
% pat1 = '*-soln.mat';
% pat2 = '*-soln.dat';
% S = struct();
% TMP1 = dir(fullfile(folder,'**\',pat1));
% N_files1 = length(TMP1);
% folders1 = cell(N_files1,1);
% 
% if (nargin==1)
%     read_again=false;
% end
% 
% if ~read_again
%     if ~isempty(TMP1)
%         % get the corresponding folders and soln (.mat) data
%         for i = 1:N_files1
%             folders1{i} = replace(TMP1(i).folder,{folder,pat1},{'',''});
%             S(i).N  = cellfun(@str2double,regexp(folders1{i},'\d*','match'));
%             S(i).N  = sqrt(prod(S(i).N));
%             file_name = fullfile(TMP1(i).folder,TMP1(i).name);
%             load(file_name,'DATA');
%             S(i).DATA = DATA;
%         end
%     end
% end
% 
% 
% TMP2 = dir(fullfile(folder,'**\',pat2));
% N_files2 = length(TMP2);
% folders2 = cell(N_files2,1);
% % get the corresponding folders and soln (.dat) data (if not already loaded)
% for i = 1:N_files2
%     folders2{i} = replace(TMP2(i).folder,{folder,pat2},{'',''});
%     % don't read in the same data again
%     if any(cellfun(@(s2)strcmp(folders2{i},s2),folders1))
%         continue
%     end
%     S(i).N  = cellfun(@str2double,regexp(folders2{i},'\d*','match'));
%     S(i).N  = sqrt(prod(S(i).N));
%     file_name = fullfile(TMP2(i).folder,TMP2(i).name);
%     % DATA = tec2mat_structured_condensed(file_name);
%     DATA = get_SENSEI_tecplot_soln_data(file_name);
%     file_name2 = replace(file_name, replace(pat2,'*',''), replace(pat1,'*','') );
%     save(file_name2,'DATA');
%     S(i).DATA = DATA;
% end
% 
% end