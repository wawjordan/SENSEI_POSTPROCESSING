%% Temp Script for calculating Richardson Extrapolation
clc; clear; close all; fclose('all');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent_dir_str = 'SENSEI_POSTPROCESSING';
path_parts = regexp(mfilename('fullpath'), filesep, 'split');
path_idx = find(cellfun(@(s1)strcmp(s1,parent_dir_str),path_parts));
parent_dir = fullfile(path_parts{1:path_idx});
addpath(genpath(parent_dir));
clear parent_dir_str path_idx path_parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;

%% for use with 'parse_and_plot_new2'
DATA_DIR='C:\Users\wajordan\Desktop\';
foldernames1 = {};
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-02_12.18.32_VAR_REC_8_1000_no_bc_IC_10_ur0-5'];% 20


foldernames1 = cellfun(@(str_b)strcat(DATA_DIR,str_b),foldernames1,UniformOutput=false);


var_list = {'RHO','U','V','P','CONS_RHO','CONS_RHO_U','CONS_RHO_V','CONS_RHO_E'};
p_hat = 2.0;
read_again = true;
use_exact  = true;

n_folders = numel(foldernames1);

for i = 1:n_folders
    fprintf('folder %d/%d:\n',i,n_folders)
    [ ~, ~, ~, ~ ] = RE_on_grid_family_data_output(foldernames1{i},p_hat,var_list,read_again,use_exact);
end