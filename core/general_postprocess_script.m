%% post-processing script
clc; clear; close all;
%%
dim   = 2;
r_fac = 2;

%% for use with 'parse_and_plot_new2'
% foldernames = {'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2024-12-12_12.31.32'};
% foldernames = [foldernames,'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2024-12-13_10.04.45'];
foldernames = {            'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_bc_constraints_IC100_Iter100_alpha_0_UR0.5_2024-12-12_12.31.32'};
foldernames = [foldernames,'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_bc_constraints_IC150_Iter150_alpha_0_UR0.5_2024-12-13_10.04.45'];
foldernames = [foldernames,'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_bc_constraints_IC_alpha_1_UR0.1_2024-12-15_12.51.01'];
foldernames = [foldernames,'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_bc_constraints_IC_alpha_1_UR0.5_2024-12-14_23.58.49'];
foldernames = [foldernames,'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_bc_constraints_IC_full_alpha_1_UR0.5_2024-12-14_11.05.56'];

foldernames = [foldernames,'C:\Users\wajordan\Desktop\SVF_subsonic_curved_bc_constraints_IC150_Iter0_2024-12-15_17.42.42'];
foldernames = [foldernames,'C:\Users\wajordan\Desktop\SVF_subsonic_curved_bc_constraints_IC150_Iter100_2024-12-15_21.24.25'];

foldernames = [foldernames,'C:\Users\wajordan\Desktop\RESULTS\AIRFOIL\CASE_3\JOUKOWSKI_C_GRID_curved_2024-12-16_12.35.28'];

foldernames = [foldernames,'C:\Users\wajordan\Desktop\Desktop_TEMP\RESULTS\AIRFOIL\CASE_3_new_investigation\JOUKOWSKI_C_GRID_curved_2025-01-22_18.23.13'];
foldernames = [foldernames,'C:\Users\wajordan\Desktop\Desktop_TEMP\RESULTS\AIRFOIL\CASE_3_new_investigation\JOUKOWSKI_C_GRID_curved_2_2025-01-22_18.23.17'];



%% effect of including boundary cells on norm
% foldernames = foldernames([1]);
% var_select    = [ 5, 6, 6 ];
% var_mask      = {[ 1, 1, 1, 1 ]};
% norm_select   = [1];
% iter_select   = {[],[100],[100]};
% layer_select  = {[],[1],[3]};
% line_fmt      = { '-', '--', ':' };
% color_spec    = {lines(4)};
% tag_fmt       = { '', '', ' (-1 layer boundary cells)' };
% legend_flag   = false;

foldernames = foldernames([9,9,10]);
var_select    = [ 3, 4, 4 ];
var_mask      = {[ 1, 0, 0, 0 ]};
norm_select   = [3];
iter_select   = {[],[0:150],[0:150]};
layer_select  = {[]};
line_fmt      = { '-', '--', '-.' };
% color_spec    = {turbo(1),jet(1),hsv(1),lines(1)};
color_spec    = {lines(4)};
tag_fmt       = { '' };
legend_flag   = false;



% foldernames = foldernames([6,6,6,7]);
% var_select    = [ 3, 4, 4, 4 ];
% var_mask      = {[ 0, 0, 1, 0 ]};
% norm_select   = [1];
% iter_select   = {[],[0],[150],[20]};
% layer_select  = {[]};
% line_fmt      = { '-', '--', '-.', ':' };
% color_spec    = {turbo(1),jet(1),hsv(1),lines(1)};
% % color_spec    = {lines(4)};
% tag_fmt       = { '' };
% legend_flag   = false;

%% effect of iterative corrections
% foldernames = foldernames([1,1,1,2,2]);
% var_select    = [ 3, 4, 4, 4, 4 ];
% var_mask      = {[ 0, 0, 0, 1 ]};
% norm_select   = [1];
% iter_select   = {[],[0],[100],[0],[150]};
% layer_select  = {[]};
% line_fmt      = { '-', '--', ':', '--', ':' };
% color_spec    = {turbo(4),jet(4),jet(4),lines(4),lines(4)};
% % color_spec    = {lines(4)};
% tag_fmt       = { '' };
% legend_flag   = false;

%% effect of iterative corrections
% foldernames = foldernames([2]);
% var_select    = [ 3, 4, 4, 4 ];
% var_mask      = {[ 1, 0, 0, 0 ]};
% norm_select   = [2];
% iter_select   = {[],[0],[1:99],[100]};
% layer_select  = {[]};
% line_fmt      = { '-', '--', '-.', ':' };
% color_spec    = {turbo(4),jet(4),lines(4),hsv(4)};
% tag_fmt       = { '' };
% legend_flag   = false;



print_ERR=false;
print_OOA=false;
err_file = 'ERR.png';
ooa_file = 'OOA.png';


%%
target_folder = fullfile('E:\wajordan\Desktop_Overflow/',string(datetime(datetime(),'Format','yyyy-MM-dd_HH.mm.ss')));


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
% set(hfig1.Children(4),'Ylim',[1e-17,1e-12])
% set(hfig1.Children(4),'Ylim',[1e-12,1e0])
% yticks(hfig1.Children(4),10.^(-12:2:0))
set(hfig1.Children(4),'Xlim',[10,10000])
% set(hfig1.Children(3),'Visible','off');
set(hfig1.Children(3),'Location','southwest');

set(hfig1.Children(2),'Ylim',[0,6])
% yticks(hfig1.Children(2),0:6)
set(hfig1.Children(2),'Xlim',[10,10000])
set(hfig1.Children(1),'Visible','off');
% set(hfig1.Children(1),'Location','southwest');


% set(hfig1.Children(3),'XlimMode','auto')
% set(hfig1.Children(2),'XLimMode','auto')
% hfig1.Children(2).Legend.Location = 'southwest';

if (print_ERR)
    exportgraphics(hfig1.Children(4),fullfile(target_folder,err_file),'Resolution',600)
end
if (print_OOA)
    exportgraphics(hfig1.Children(2),fullfile(target_folder,ooa_file),'Resolution',600)
end

