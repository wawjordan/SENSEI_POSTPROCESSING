%% Parsing SVF-airfoil data (12/15/2025)
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

dim   = 2;
r_fac = 2;

%% for use with 'parse_and_plot_new2'
DATA_DIR='C:\Users\wajordan\Desktop\';
foldernames1 = {};
foldernames1 = [foldernames1,'SVF_2025-12-15_15.31.44_KEXACT_EXHO_EVAL=T_REC=T'];% 1
foldernames1 = [foldernames1,'SVF_2025-12-15_15.37.24_VAR_REC_4_20_EXHO_EVAL=T_REC=T'];% 2
foldernames1 = [foldernames1,'SVF_2025-12-15_15.42.54_VAR_REC_4_100_EXHO_EVAL=T_REC=T'];% 3
foldernames1 = [foldernames1,'SVF_2025-12-15_15.47.26_VAR_REC_4_1000_EXHO_EVAL=T_REC=T'];% 4
foldernames1 = [foldernames1,'SVF_2025-12-15_16.03.52_VAR_REC_4_2000_EXHO_EVAL=T_REC=T'];% 5
foldernames1 = [foldernames1,'SVF_2025-12-15_16.13.29_VAR_REC_4_4000_EXHO_EVAL=T_REC=T'];% 6
foldernames1 = [foldernames1,'SVF_2025-12-15_16.19.42_VAR_REC_8_1000_EXHO_EVAL=T_REC=T'];% 7
% foldernames1 = [foldernames1,'SVF_2025-12-15_16.27.24_VAR_REC_8_1000_EXHO_EVAL=T_REC=T_OLDER'];% 8
% foldernames1 = [foldernames1,'SVF_2025-12-15_16.30.59_VAR_REC_8_1000_EXHO_EVAL=T_REC=T_OLDEST'];% 9
foldernames1 = [foldernames1,'SVF_2025-12-15_17.20.13_VAR_REC_8_1000_EXHO_EVAL=T_REC=T_2575a92'];% 8
foldernames1 = [foldernames1,'SVF_2025-12-15_17.25.43_VAR_REC_1_10_EXHO_EVAL=T_REC=T_2575a92'];% 9
foldernames1 = [foldernames1,'SVF_2025-12-15_17.28.37_VAR_REC_ALL_10_EXHO_EVAL=T_REC=T_2575a92'];% 10
foldernames1 = [foldernames1,'SVF_2025-12-15_17.59.59_VAR_REC_ALL_2000_EXHO_EVAL=T_REC=T_2575a92'];% 11
foldernames1 = [foldernames1,'SVF_2025-12-15_18.11.47_VAR_REC_4_10_EXHO_EVAL=T_REC=T_200IC_2575a92'];% 12
foldernames1 = [foldernames1,'SVF_2025-12-15_18.18.13_KEXACT_EXHO_EVAL=T_REC=T_200IC_2575a92'];% 13
foldernames1 = [foldernames1,'SVF_2025-12-15_18.28.16_VAR_REC_1_1_EXHO_EVAL=T_REC=T_200IC_2575a92'];% 14
foldernames1 = [foldernames1,'SVF_2025-12-15_18.35.44_VAR_REC_ALL_1_EXHO_EVAL=T_REC=T_200IC_2575a92'];% 15
foldernames1 = [foldernames1,'SVF_2025-12-15_18.55.10_VAR_REC_1_1000_EXHO_EVAL=T_REC=T_200IC_2575a92'];% 16
foldernames1 = [foldernames1,'SVF_2025-12-16_12.19.22_VAR_REC_1_1000_EXHO_EVAL=T_REC=T_200IC_yes_bc_2575a92'];% 17
foldernames1 = [foldernames1,'SVF_2025-12-18_11.36.46_VAR_REC_1_1000_EXHO_EVAL=T_REC=T_200IC_yes_bc2_2575a92'];% 18
foldernames1 = [foldernames1,'SVF_2025-12-18_11.51.30_KEXACT_EXHO_EVAL=T_REC=T_200IC_yes_bc_2575a92'];% 19
foldernames1 = [foldernames1,'SVF_2025-12-18_13.49.39_KEXACT_EXHO_EVAL=T_REC=T_200IC_yes_bc_geo4_2575a92'];% 20
foldernames1 = [foldernames1,'SVF_2025-12-18_14.06.56_KEXACT3_EXHO_EVAL=T_REC=T_200IC_yes_bc_2575a92'];% 21


foldernames1 = cellfun(@(str_b)strcat(DATA_DIR,str_b),foldernames1,UniformOutput=false);

foldernames = foldernames1([13,19,21]);
var_select    = [ 4, 4, 4 ];
var_mask      = {[ 1, 1, 1, 1 ]};
norm_select   = [3];
iter_select   = {[]};
layer_select  = {[]};
line_fmt      = { '-', ':', '--' };
color_spec    = {lines(4)};
tag_fmt       = { '' };
legend_flag   = true;
post_plot_commands = {"set(hfig1.Children(4),'Ylim',[1e-11,1e-3]);",...
                      "yticks(hfig1.Children(4),10.^(-11:1:-3));",  ...
                      "set(hfig1.Children(4),'Xlim',[10,1000])",    ...
                      "set(hfig1.Children(3),'Location','southwest');",...
                      "set(hfig1.Children(2),'Ylim',[0,5]);",...
                      "set(hfig1.Children(2),'Xlim',[10,1000])",...
                      "set(hfig1.Children(1),'Visible','off');"};
% post_plot_commands = {"set(hfig1.Children(4),'Ylim',[1e-5,1e0]);",...
%                       "yticks(hfig1.Children(4),10.^(-5:1:0));",  ...
%                       "set(hfig1.Children(4),'Xlim',[10,1000])",    ...
%                       "set(hfig1.Children(3),'Location','southwest');",...
%                       "set(hfig1.Children(2),'Ylim',[0,5]);",...
%                       "set(hfig1.Children(2),'Xlim',[10,1000])",...
%                       "set(hfig1.Children(1),'Visible','off');"};
print_ERR=false;
print_OOA=false;
target_folder = 'C:\Users\wajordan\Desktop\';
err_file = '';
ooa_file = '';

[hfig1,DE_test] = parse_and_plot_new2(dim,r_fac, foldernames,          ...
                                                          var_select,   ...
                                                          var_mask,     ...
                                                          norm_select,  ...
                                                          iter_select,  ...
                                                          layer_select, ...
                                                          tag_fmt,      ...
                                                          line_fmt,     ...
                                                          color_spec,   ...
                                                          legend_flag );
cellfun(@eval,post_plot_commands)

if (print_ERR)
    exportgraphics(hfig1.Children(4),fullfile(target_folder,err_file),'Resolution',600)
end
if (print_OOA)
    exportgraphics(hfig1.Children(2),fullfile(target_folder,ooa_file),'Resolution',600)
end