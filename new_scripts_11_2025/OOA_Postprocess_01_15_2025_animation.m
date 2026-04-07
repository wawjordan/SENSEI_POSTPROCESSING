%% Parsing KT-airfoil data (11/23/2025)
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
DATA_DIR='C:\Users\wajordan\Desktop\CASES\';
foldernames1 = {             'ALPHA_0_JOUKOWSKI_C_GRID_curved_2025-11-23_13.28.07_K_EXACT'};% 1
foldernames1 = [foldernames1,'ALPHA_0_JOUKOWSKI_C_GRID_curved_2025-11-23_12.47.49_VAR_REC_8_100'];% 2
foldernames1 = [foldernames1,'ALPHA_0_JOUKOWSKI_C_GRID_curved_2025-11-23_14.30.04_VAR_REC_8_1000'];% 3
foldernames1 = [foldernames1,'ALPHA_0_JOUKOWSKI_C_GRID_curved_2025-11-23_15.35.53_VAR_REC_ALL_100']; % 4
foldernames1 = [foldernames1,'ALPHA_0_JOUKOWSKI_C_GRID_curved_2025-11-23_15.03.36_VAR_REC_ALL_1000'];% 5
foldernames1 = [foldernames1,'ALPHA_0_JOUKOWSKI_C_GRID_curved_2025-11-23_15.57.49_VAR_REC_ALL_100_USE_HO_GRID_RECONSTRUCT=T'];% 6
foldernames1 = [foldernames1,'ALPHA_0_JOUKOWSKI_C_GRID_curved_2025-11-23_16.13.48_VAR_REC_ALL_100_no_bc'];% 7
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-11-26_15.23.56_K_EXACT_no_bc'];% 8
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-11-26_15.40.39_K_EXACT'];% 9
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-11-23_16.46.56_VAR_REC_ALL_100_no_bc'];% 10
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-11-26_16.20.39_VAR_REC_ALL_1000_no_bc'];% 11
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-11-26_15.54.04_VAR_REC_ALL_100_bc'];% 12

foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-04_12.39.32_KEXACT_yes_bc_IC_10_ur0-5'];% 13


foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-01_12.50.36_VAR_REC_ALL_1000_USE_HO_GRID_RECONSTRUCT=T_no_bc_IC_10_no_ur'];% 14
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-01_23.04.15_VAR_REC_ALL_1000_no_bc_IC_10_ur0-5'];% 15
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-02_10.32.10_VAR_REC_ALL_10_no_bc_IC_10_ur0-5'];% 16

foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-02_19.07.11_VAR_REC_1_1000_no_bc_IC_10_ur0-5'];% 17
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-02_16.26.38_VAR_REC_2_1000_no_bc_IC_10_ur0-5'];% 18
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-02_14.01.27_VAR_REC_4_1000_no_bc_IC_10_ur0-5'];% 19
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-02_12.18.32_VAR_REC_8_1000_no_bc_IC_10_ur0-5'];% 20

foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-03_20.39.41_VAR_REC_8_10_yes_bc_IC_10_ur0-5'];% 21
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-04_11.01.38_VAR_REC_8_100_yes_bc_IC_10_ur0-5'];% 22
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-03_11.55.12_VAR_REC_8_1000_yes_bc_IC_10_ur0-5'];% 23
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-05_11.19.30_VAR_REC_1_100_yes_bc_IC_10_ur0-5'];% 24
% abe23e21f0a6d7f1466407415ad2ebb6f5925edc


foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-02_17.52.48_VAR_REC_8_100_no_bc_IC_200_ur0-5'];% 25
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-08_13.39.39_VAR_REC_ALL_1000_no_bc_IC_200_ur0-5'];% 26


%% New compile (12d0c99f7af602e3eac64ad5bde538f3262efefb)
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-05_13.20.34_VAR_REC_1_100_yes_bc_IC_10_ur0-5_NEW'];% 27
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-05_18.06.19_VAR_REC_1_100_yes_bc_IC_10_ur0-5_NEW_10_iter'];% 28
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-08_13.32.09_VAR_REC_2_100_yes_bc_IC_10_ur0-5_NEW_10_iter'];% 29
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-09_13.43.19_ORDER_3_VAR_REC_1_100_yes_bc_IC_10_ur0-5_NEW_10_iter'];% 30
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-09_19.44.01_ORDER_3_VAR_REC_1_100_no_bc_IC_10_ur0-5_NEW_10_iter'];% 31
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_newer_2025-12-09_21.41.18_ORDER_3_VAR_REC_1_100_no_bc_IC_10_ur0-5_NEWER_10_iter'];% 32
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-10_10.50.20_ORDER_3_VAR_REC_1_100_no_bc_IC_10_ur0-5_NEW_10_iter_GEO4'];% 33

foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2026-01-15_02.00.10_CURVED_P4_q6_q6_hoex_TT_BC_VARLAY2_ITER2000_IC10'];% 34
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2026-01-15_03.10.03_CURVED_P4_q6_q6_hoex_TT_BC_VARLAY8_ITER2000_IC10'];% 35

foldernames1 = cellfun(@(str_b)strcat(DATA_DIR,str_b),foldernames1,UniformOutput=false);

%% variational (8 layers) vs kexact w/constraints
foldernames = foldernames1([13,13,35]); tag_fmt = { '', '(Old Rec.)', '(New Rec.)' };
var_select    = [ 3, 4, 4 ];
var_mask      = {[ 1, 1, 1, 1 ]};
norm_select   = [1];
layer_select  = {[]};
line_fmt      = { '-', ':', '--' };
color_spec    = {lines(4)};
legend_flag   = true;

post_plot_commands = {"set(hfig1.Children(4),'Ylim',[1e-11,1e-3]);",...
                      "yticks(hfig1.Children(4),10.^(-11:1:-3));",  ...
                      "set(hfig1.Children(4),'Xlim',[10,1000])",    ...
                      "set(hfig1.Children(3),'Location','southwest');",...
                      "set(hfig1.Children(2),'Ylim',[0,5]);",...
                      "set(hfig1.Children(2),'Xlim',[10,1000])",...
                      "set(hfig1.Children(1),'Visible','off');"};
print_ERR=true;
print_OOA=true;
target_folder = 'C:\Users\wajordan\Desktop\SciTech_Plots\ERR_norms_airfoil';
for j = 0:10
iter_select = {[j],[j],[j]};
err_file = sprintf('ERR_%0.2d.png',j);
ooa_file = sprintf('OOA_%0.2d.png',j);

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
label_iteration(hfig1.Children(4),j);
label_iteration(hfig1.Children(2),j);
if (print_ERR)
    exportgraphics(hfig1.Children(4),fullfile(target_folder,err_file),'Resolution',600)
end
if (print_OOA)
    exportgraphics(hfig1.Children(2),fullfile(target_folder,ooa_file),'Resolution',600)
end

end

function label_iteration(ax,iter)

xpos = 0.4;
ypos = 0.97;
txt = sprintf('Iterative Correction: %d',iter);

text(ax,xpos,ypos,txt,'Interpreter','latex','Units','Normalized',...
                      'VerticalAlignment','top')
end