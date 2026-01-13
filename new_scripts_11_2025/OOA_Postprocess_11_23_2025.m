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


foldernames1 = cellfun(@(str_b)strcat(DATA_DIR,str_b),foldernames1,UniformOutput=false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% variational (8 layers) vs kexact w/o/constraints
% foldernames = foldernames1([20,9,20]); iter_select = {[],[10],[10]}; tag_fmt = { '', '(k-exact)', '(variational [8], 1000 iter)' };

%% variational (all layers) vs kexact w/o/constraints
% foldernames = foldernames1([15,9,15]); iter_select = {[],[10],[10]}; tag_fmt = { '', '(k-exact)', '(variational [all], 1000 iter)' };


%% variational (1 v 8 layers) w/constraints
% foldernames = foldernames1([24,24,22]); iter_select = {[],[10],[10]}; tag_fmt = { '', '(variational [1])', '(variational [8])' };



% foldernames = foldernames1([13,13,24]); iter_select = {[],[0],[0]}; tag_fmt = { '', '(k-exact)', '(variational [1])' };

% foldernames = foldernames1([17,17,24]); iter_select = {[],[10],[10]}; tag_fmt = { '', '(no constraints)', '(yes constraints)' };

% foldernames = foldernames1([17,8,17]); iter_select = {[],[0],[0]}; tag_fmt = { '', '(k-exact)', '(variational [1])' };


% foldernames = foldernames1([24,24,27]); iter_select = {[],[0:10],[0:10]}; tag_fmt = { '', '(OLD)', '(NEW)' };

%% variational (1 v 2 layers) New compile
% foldernames = foldernames1([28,28,29]); iter_select = {[],[1],[1]}; tag_fmt = { '', '(1)', '(2)' };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% variational (8 layers) vs kexact w/constraints
foldernames = foldernames1([13,13,20]); iter_select = {[],[],[]}; tag_fmt = { '', '(Old Rec.)', '(New Rec.)' };

%% variational (all layers) vs kexact w/constraints
% foldernames = foldernames1([26,26,26]); iter_select = {[],[0],[1:5:200]}; tag_fmt = { '', '(variational [all], 1000 iter (10))', '(variational [all], 1000 iter (200))' };


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% variational (8 layers) vs variational (all) w/constraints [200 IC]
% foldernames = foldernames1([25,25,25]); iter_select = {[],[0],[10]}; tag_fmt = { '', '(variational [8], 100 iter)', '(variational [8], 100 iter, 200 IC)' };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% variational (8 layers w/o constraints) vs variational (8 w/constraints)
% foldernames = foldernames1([20,20,23]); iter_select = {[],[0:10],[0:10]}; tag_fmt = { '', '(variational [8], no constraints)', '(variational [8], constraints)' };

% foldernames = foldernames1([30,30,28]); iter_select = {[],[],[0:10]}; tag_fmt = { '', '(r=3)', '(r=4)' };
% foldernames = foldernames1([31,31,26]); iter_select = {[],[],[0:200]}; tag_fmt = { '', '(r=3)', '(r=4)' };
% foldernames = foldernames1([30,31,30]); iter_select = {[],[0:2],[0:2]}; tag_fmt = { '', '(r=3)', '(r=3 w/constraints)' };

% foldernames = foldernames1([31,31,32]); iter_select = {[],[0],[0]}; tag_fmt = { '', '(OLD)', '(NEW)' };
% foldernames = foldernames1([32,32,33]); iter_select = {[],[0],[0]}; tag_fmt = { '', '(2)', '(4)' };
% 
var_select    = [ 3, 4, 4 ];
var_mask      = {[ 0, 1, 0, 1 ]};
norm_select   = [3];
layer_select  = {[]};
line_fmt      = { '-', ':', '--' };
color_spec    = {lines(4)};
legend_flag   = true;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Effect of num iterations on var rec (8 layers)
% foldernames = foldernames1([21,21,22,23]);
% var_select    = [ 3, 4, 4, 4 ];
% var_mask      = {[ 1, 1, 1, 1 ]};
% norm_select   = [2];
% iter_select   = {[],[],[],[]};
% layer_select  = {[]};
% line_fmt      = { '-', '--', '-.', ':' };
% % color_spec    = {turbo(1),jet(1),hsv(1),lines(1)};
% color_spec    = {lines(4)};
% tag_fmt       = { '', '(variational [8], 10 iter)', '(variational [8], 100 iter)', '(variational [8], 1000 iter)' };
% legend_flag   = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Effect of num layers on var rec (1000 iter)
% foldernames = foldernames1([17,18,19,20]);
% var_select    = [ 4, 4, 4, 4 ];
% var_mask      = {[ 1, 1, 1, 1 ]};
% norm_select   = [3];
% iter_select   = {[],[],[],[]};
% layer_select  = {[]};
% line_fmt      = { '-', '--', '-.', ':' };
% color_spec    = {lines(4)};
% tag_fmt       = { '(variational [1])', '(variational [2])', '(variational [4])', '(variational [8])' };
% legend_flag   = true;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% foldernames = foldernames1([1,2]);
% var_select    = [ 6, 6 ];
% var_mask      = {[ 1, 1, 1, 1 ]};
% norm_select   = [2];
% iter_select   = {[]};
% layer_select  = {[]};
% line_fmt      = { '-', '--' };
% color_spec    = {lines(4)};
% tag_fmt       = { '' };
% legend_flag   = true;
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
print_ERR=true;
print_OOA=true;
target_folder = 'C:\Users\wajordan\Desktop\CCAS_Annual_Review_Plots\OOA\L1_NORM';
err_file = 'ERR_L1_U_and_P_only_0_IC.png';
ooa_file = 'OOA_L1_U_and_P_only_0_IC.png';

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

% set(hfig1.Children(4),'Ylim',[1e-17,1e-12])
% set(hfig1.Children(4),'Ylim',[1e-12,1e0])
% yticks(hfig1.Children(4),10.^(-12:2:0))
% set(hfig1.Children(4),'Xlim',[10,10000])
% set(hfig1.Children(3),'Visible','off');
% set(hfig1.Children(3),'Location','southwest');
% set(hfig1.Children(2),'Ylim',[0,6])
% yticks(hfig1.Children(2),0:6)
% set(hfig1.Children(2),'Xlim',[10,10000])
% set(hfig1.Children(1),'Visible','off');
% set(hfig1.Children(1),'Location','southwest');
% set(hfig1.Children(3),'XlimMode','auto')
% set(hfig1.Children(2),'XLimMode','auto')
% hfig1.Children(2).Legend.Location = 'southwest';