%% Parsing KT-airfoil data (04/02/2026)
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
foldernames1 = {};
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2026-01-15_02.00.10_CURVED_P4_q6_q6_hoex_TT_BC_VARLAY2_ITER2000_IC10'];% 1
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2026-01-15_03.10.03_CURVED_P4_q6_q6_hoex_TT_BC_VARLAY8_ITER2000_IC10'];% 2
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2026-02-13_12.01.06_ORDER_4_bc_IC_10_100_iter_GEO4'];% 3
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_limited_2026-02-13_12.00.09_ORDER_4_bc_IC_10_100_iter_GEO4'];% 4

foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_04_01_2026-04-02_12.20.09'];        % 5
foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_04_01_regress_2026-04-02_12.21.04'];% 6

foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_clustered_04_02_regress_2026-04-02_14.54.19'];% 7
foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_clustered_04_02_2026-04-02_14.54.15'];% 8

foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_clustered_04_02_regress_2026-04-02_19.10.44'];% 9
foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_clustered_04_02_2026-04-02_19.10.39'];% 10

foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_clustered_04_02_regress_2026-04-03_11.01.28'];% 11
foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_clustered_04_02_2026-04-03_10.47.17'];% 12


foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_v2_04_08_regress_2026-04-09_11.38.29']; % 13
foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_v2_04_08_2026-04-09_11.38.25'        ]; % 14

foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_v3_04_09_regress_2026-04-09_18.12.43'        ]; % 15
foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_v3_04_09_2026-04-09_18.12.39'        ]; % 16

foldernames1 = [foldernames1,  'JOUKOWSKI_C_GRID_curved_v3_04_09_regress_2026-04-10_11.40.43'      ]; % 17
foldernames1 = [foldernames1,  'JOUKOWSKI_C_GRID_curved_v3_04_09_2026-04-10_11.40.38'      ]; % 18

foldernames1 = [foldernames1,  'REGRESS_4_JOUKOWSKI_C_GRID_curved_v3_04_09_regress_2026-05-13_12.36.01'      ]; % 19

foldernames1 = [foldernames1,  'ITER_0_JOUKOWSKI_C_GRID_curved_v3_04_09_regress_2026-05-14_20.02.21'      ]; % 20






foldernames1 = cellfun(@(str_b)strcat(DATA_DIR,str_b),foldernames1,UniformOutput=false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% foldernames = foldernames1([44,44,45]); iter_select = {[],[150],[150]}; tag_fmt = { '', '(k-exact)', '(cweno)' };
% var_select    = [ 3, 4, 4 ];
% var_mask      = {[ 1, 1, 1, 1 ]};
% norm_select   = [1];
% layer_select  = {[]};
% line_fmt      = { '-', '-.', '--' };
% color_spec    = {lines(4)};
% legend_flag   = true;

% foldernames = foldernames1([14,13,14]); iter_select = {[],[],[]}; tag_fmt = { '', '(k-exact)', '(cweno)' };
% var_select    = [ 3, 4, 4 ];
% var_mask      = {[ 1, 1, 1, 1 ]};
% norm_select   = [1];
% layer_select  = {[]};
% line_fmt      = { '-', '-.', '--' };
% color_spec    = {lines(4)};
% legend_flag   = true;
% foldernames = foldernames1([18,17,18]); iter_select = {[],[],[]}; tag_fmt = { '', '(k-exact)', '(cweno)' };
% var_select    = [ 5, 6, 6 ];
% var_mask      = {[ 1, 1, 1, 1 ]};
% norm_select   = [3];
% layer_select  = {[]};
% line_fmt      = { '-', '-.', '--' };
% color_spec    = {lines(4)};
% legend_flag   = true;



% foldernames = foldernames1([17,19,17,19]); iter_select = {[],[],[],[]}; tag_fmt = { '(unlimited)', '(limited)', '(unlimited)', '(limited)' };
% var_select    = [ 3, 3, 4, 4 ];
% var_mask      = {[ 1, 1, 1, 1 ]};
% norm_select   = [1];
% layer_select  = {[]};
% line_fmt      = { ':','-', '-.', '--' };
% color_spec    = {lines(4)};
% legend_flag   = true;

% foldernames = foldernames1([17,19]); iter_select = {[],[]}; tag_fmt = { '' };
% var_select    = [ 4, 4 ];
% var_mask      = {[ 1, 1, 1, 1 ]};
% norm_select   = [1];
% layer_select  = {[]};
% line_fmt      = { '-', '--' };
% color_spec    = {lines(4)};
% legend_flag   = true;

foldernames = foldernames1([19,19,20]); iter_select = {[],[10],[10]}; tag_fmt = { '' };
var_select    = [ 3, 4, 4 ];
var_mask      = {[ 1, 1, 1, 1 ]};
norm_select   = [1];
layer_select  = {[]};
line_fmt      = { '-', '--', ':' };
color_spec    = {lines(4)};
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