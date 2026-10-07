%% Parsing KT-airfoil data (09/17/2026)
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


% foldernames1 = [foldernames1,'ms1_test_09172026_2026-09-17_15.18.34']; % bc=mms, nskip=4, 1 layer curv,   0 iter, UR=0.1, nl_iter=0,  bc_con=T, cbc=T
% foldernames1 = [foldernames1,'ms1_test_09172026_2026-09-17_17.33.51']; % bc=mms, nskip=4, 1 layer curv,  50 iter, UR=0.1, nl_iter=50, bc_con=T, cbc=T
% foldernames1 = [foldernames1,'ms1_test_09172026_2026-09-20_17.23.59']; % bc=mms, nskip=4, 1 layer curv,  50 iter, UR=0.1, nl_iter=50, bc_con=F, cbc=T
foldernames1 = [foldernames1,'ms1_test_09172026_2026-09-21_07.23.16']; % bc=mms, nskip=4, 1 layer curv, 200 iter, UR=0.5, nl_iter=50, bc_con=T, cbc=T

foldernames1 = [foldernames1,'MOD_ms1_test_09232026_mms_2026-09-24_00.55.36']; % bc=mms, nskip=4, 1 layer curv, 200 iter, UR=0.5, nl_iter=50, bc_con=T, cbc=T
foldernames1 = [foldernames1,'MOD2_ms1_test_09232026_mms_2026-09-24_10.52.21']; % bc=mms, nskip=4, 1 layer curv, 200 iter, UR=0.5, nl_iter=50, bc_con=T, cbc=T

foldernames1 = [foldernames1,'NO_CBC_ms1_test_09232026_mms_2026-09-24_14.38.02']; % bc=mms, nskip=4, 1 layer curv, 200 iter, UR=0.5, nl_iter=50, bc_con=T, cbc=F

% foldernames1 = [foldernames1,'ms1_test_09212026_2026-09-21_10.27.17']; % bc=slip, nskip=4,-1 layer curv, 200 iter, UR=0.5, nl_iter=50, bc_con=T, cbc=T
% foldernames1 = [foldernames1,'ms1_test_09212026_2026-09-21_11.03.39']; % bc=slip, nskip=4,-1 layer curv, 200 iter, UR=0.5, nl_iter=50, bc_con=T, cbc=T, LETE
% foldernames1 = [foldernames1,'ms1_test_09212026_2026-09-22_10.48.13']; % bc=slip, nskip=4,-1 layer curv, 500 iter, UR=0.5, nl_iter=50, bc_con=T, cbc=F
% foldernames1 = [foldernames1,'ms1_test_09212026_2026-09-22_14.53.26']; % bc=slip, nskip=2,-1 layer curv, 200 iter, UR=0.5, nl_iter=50, bc_con=T, cbc=T
% foldernames1 = [foldernames1,'ms1_test_09212026_2026-09-23_01.49.15']; % bc=slip, nskip=4,-1 layer curv, 200 iter, UR=0.5, nl_iter=50, bc_con=T, cbc=T, REC=3
% foldernames1 = [foldernames1,'ms1_test_09212026_2026-09-23_15.37.13']; % bc=slip, nskip=4,-1 layer curv, 200 iter, UR=0.5, nl_iter=50, bc_con=T, cbc=T, REC=4
% foldernames1 = [foldernames1,'ms1_test_09212026_2026-09-23_09.11.33']; % bc=slip, nskip=4,-1 layer curv, 200 iter, UR=0.5, nl_iter=50, bc_con=T, cbc=T, REC=5


foldernames1 = [foldernames1,'stag_flow_test_10_01_2026_2026-10-02_18.19.15'];

foldernames1 = [foldernames1,'stag2_test_10_05_2026_2026-10-05_10.12.37'];

foldernames1 = [foldernames1,'parabola_test_AR_1_2026-10-07_18.40.06'];


foldernames1 = cellfun(@(str_b)strcat(DATA_DIR,str_b),foldernames1,UniformOutput=false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% foldernames = foldernames1([1,2,4]); iter_select = {[0],[0],[]}; tag_fmt = { '' };
% var_select    = [ 4, 4, 4 ];
% var_mask      = {[ 1, 1, 1, 1 ]};
% norm_select   = [1];
% layer_select  = {[]};
% line_fmt      = { '-', '--',':'};
% color_spec    = {lines(4),lines(4),hsv(4)};
% legend_flag   = true;

foldernames = foldernames1([7,7]); iter_select = {[0],[]}; tag_fmt = { '' };
var_select    = [ 3, 4 ];
var_mask      = {[ 1, 1, 1, 1 ]};
norm_select   = [1];
layer_select  = {[]};
line_fmt      = { '-', '--'};
color_spec    = {lines(4)};
legend_flag   = true;


post_plot_commands = {"set(hfig1.Children(4),'Ylim',[1e-11,1e-3]);",...
                      "yticks(hfig1.Children(4),10.^(-11:1:-3));",  ...
                      "set(hfig1.Children(4),'Xlim',[10,1000])",    ...
                      "set(hfig1.Children(3),'Location','southwest');",...
                      "set(hfig1.Children(2),'Ylim',[0,7]);",...
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