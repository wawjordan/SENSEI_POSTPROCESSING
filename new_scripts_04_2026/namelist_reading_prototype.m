%% SENSEI Namelist Parsing Prototype (05/15/2026)
clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent_dir_str = 'SENSEI_POSTPROCESSING';
path_parts = regexp(mfilename('fullpath'), filesep, 'split');
path_idx = find(cellfun(@(s1)strcmp(s1,parent_dir_str),path_parts));
parent_dir = fullfile(path_parts{1:path_idx});
addpath(genpath(parent_dir));
clear parent_dir_str path_idx path_parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;


file = 'C:\Users\wajordan\Desktop\CASES\ITER_0_JOUKOWSKI_C_GRID_curved_v3_04_09_regress_2026-05-14_20.02.21\_kt4097x1025\sensei.nml';

tmp = extractBetween(fileread(file),'GRID_NSKIP',newline);

file_text = fileread(file);
out_text = strip_comments(file_text);
S = get_sections(out_text);
% strip comments

extracted = extractBetween(file_text, '&', '/');
extracted = extractBetween(file_text, '&', '/');

extractBetween(extracted{6},newline,',');
extractBefore(extracted{6},newline);
% function read_nml_section()

function out_text = strip_comments(in_text)
tmp1 = extractBetween(in_text, '!', newline,'Boundaries','inclusive');
tmp1 = cellfun(@(tmp1)strrep(tmp1,newline,''),tmp1,'UniformOutput',false);
out=cellfun(@(pat)strrep(in_text,pat,''),tmp1,'UniformOutput',false);
out_text = [out{:}];
end

function S = get_sections(in_text)
tmp = extractBetween(in_text, '&', '/');
S = struct();
for i = 1:numel(tmp)
    header = extractBefore(tmp{i},newline);
    entries = extractBetween(tmp{i},newline,',');
    % S.(header) = entries;
    S.(header) = struct();
    tmp_entry_names = cell(numel(entries),1);
    tmp_entries     = cell(numel(entries),1);
    mask            = false(numel(entries),1);
    for j = 1:numel(entries)
        tmp_entry_names{j} = extractBefore(entries{j},'=');
        tmp_entry_names{j} = strrep(tmp_entry_names{j},' ','');
        tmp_entry{j} = extractAfter(entries{j},'=');
        tmp_entry{j} = strrep(tmp_entry{j},' ','');
        mask(j) = ~contains(entries{j},'=');
    end

        entry_name = extractBefore(entries{j},'=');
    %     entry_name = strrep(entry_name,' ','');
    %     entry = extractAfter(entries{j},'=');
    %     entry = strrep(entry,' ','');
    %     S.(header).(entry_name) = entry;
    % end
end
end

% function out_text = strip_trailing_whitespace(in_text)
% end
% function [DATA] = read_SENSEI_nml_file(file_name)
% fid = fopen(file_name,'r');
% 
% % the next line should be the (grid size, 1D grid size, time)
% line = fgetl(fid);
% tmp = sscanf(line,'%f');
% N = tmp(1);
% 
% if (length(tmp)==3)
%     t = tmp(3);
% else
%     t = [];
% end
% 
% line = fgetl(fid);
% tmp = sscanf(line,'%f');
% N_eqns = length(tmp);
% dat = zeros(N_eqns,3);
% dat(:,1) = tmp;
% 
% for i = 2:3
%     line = fgetl(fid);
%     tmp = sscanf(line,'%f');
%     dat(:,i) = tmp;
% end
% 
% fclose(fid);
% end