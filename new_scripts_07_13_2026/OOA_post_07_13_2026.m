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
% foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_v2_04_08_regress_regress_2026-07-13_14.26.26'];
% foldernames1 = [foldernames1,'KT_AR_1_1_P2_ALPHA_05_2026-07-14_00.34.46'];
% foldernames1 = [foldernames1,'KT_AR_1_10_P2_ALPHA_05_2026-07-14_00.35.35'];
% foldernames1 = [foldernames1,'KT_AR_1_1_P4_ALPHA_05_2026-07-20_10.39.40'];
% foldernames1 = [foldernames1,'KT_AR_1_10_P4_ALPHA_05_2026-07-20_10.40.03'];
% foldernames1 = [foldernames1,'KT_AR_1_1_P4_ALPHA_00_2026-07-20_19.05.51'];
% foldernames1 = [foldernames1,'KT_AR_1_10_P4_ALPHA_00_2026-07-20_19.06.05'];
% foldernames1 = [foldernames1,'KT_AR_1_1_P4_ALPHA_05A_2026-07-20_20.45.05'];
% foldernames1 = [foldernames1,'KT_AR_1_1_P4_ALPHA_05A_N20_2026-07-20_22.23.14'];
% foldernames1 = [foldernames1,'KT_AR_1_1_TE_10_P4_ALPHA_05_2026-07-20_22.15.53'];
% foldernames1 = [foldernames1,'KT_AR_1_1_TE_10_P4_ALPHA_05_NEW_2026-07-21_13.39.19'];
% foldernames1 = [foldernames1,'KT_AR_1_1_P4_ALPHA_05A_D1_2026-07-21_16.31.38'];
% foldernames1 = [foldernames1,'KT_AR_1_1_P4_ALPHA_05A_D02_2026-07-22_14.15.01'];
% foldernames1 = [foldernames1,'KT_AR_1_1_P4_ALPHA_05A_D005_2026-07-22_14.15.18'];
% foldernames1 = [foldernames1,'KT_AR_1_1_P4_ALPHA_00_NEW_2026-07-23_10.23.07'];
% foldernames1 = [foldernames1,'KT_AR_1_1_P4_ALPHA_05_NEW_2026-07-23_20.13.17'];
foldernames1 = [foldernames1,'KT_AR_1_1_TE_10_P4_ALPHA_05_NEW_2026-07-24_11.46.58'];
foldernames1 = [foldernames1,'KT_AR_1_1_TE_10_P4_ALPHA_05_NEW_PW_2026-07-24_18.46.06'];



foldernames1 = cellfun(@(str_b)strcat(DATA_DIR,str_b),foldernames1,UniformOutput=false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% foldernames = foldernames1([1,1,2]); iter_select = {[],[30],[30]}; tag_fmt = { '' };
% var_select    = [ 3, 4, 4];
% var_mask      = {[ 1, 1, 1, 1 ]};
% norm_select   = [3];
% layer_select  = {[]};
% line_fmt      = { '-', '--', ':'};
% color_spec    = {lines(4)};
% legend_flag   = true;

foldernames = foldernames1([1,1,2]); iter_select = {[],[],[]}; tag_fmt = { '' };
var_select    = [ 3, 4, 4 ];
var_mask      = {[ 1, 1, 1, 1 ]};
norm_select   = [3];
layer_select  = {[]};
line_fmt      = { '-', '--', ':'};
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