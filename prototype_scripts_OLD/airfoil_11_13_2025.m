%% post-processing script
% kt airfoil 11/13/2025

clc; clear; close all;
%%
dim   = 2;
r_fac = 2;

%% for use with 'parse_and_plot_new2'
% foldernames1 = {             'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2_2025-11-13_12.43.22'};                 % 1
% foldernames1 = [foldernames1,'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2_2025-11-13_12.33.34_new_k_exact'];     % 2

foldernames1 = {             'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2_2025-11-13_12.06.39_old_k_exact'};          % 1
foldernames1 = [foldernames1,'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2_2025-11-13_12.33.34_new_k_exact'];          % 2 
foldernames1 = [foldernames1,'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2_2025-11-13_12.43.22_new_var_rec_2_layers']; % 3 
foldernames1 = [foldernames1,'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2_2025-11-13_12.57.27_new_var_rec_8_layers']; % 4
foldernames1 = [foldernames1,'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2_2025-11-13_13.42.58_new_var_rec_8_layers_reuse']; % 5

foldernames = foldernames1([1,2,5]);
var_select    = [ 3, 4, 4 ];
var_mask      = {[ 1, 1, 1, 1 ]};
norm_select   = [1];
iter_select   = {[],[],[]};
layer_select  = {[]};
line_fmt      = { '-', '--', ':' };
% color_spec    = {turbo(1),jet(1),hsv(1),lines(1)};
color_spec    = {lines(4)};
tag_fmt       = { '', '(1)', '(2)' };
legend_flag   = true;
post_plot_commands = {"set(hfig1.Children(4),'Ylim',[1e-10,1e-3]);",...
                      "yticks(hfig1.Children(4),10.^(-10:1:-3));",  ...
                      "set(hfig1.Children(4),'Xlim',[10,1000])",    ...
                      "set(hfig1.Children(3),'Location','southwest');",...
                      "set(hfig1.Children(2),'Ylim',[0,5]);",...
                      "set(hfig1.Children(2),'Xlim',[10,1000])",...
                      "set(hfig1.Children(1),'Visible','off');"};
print_ERR=false;
print_OOA=false;


target_folder = 'C:\Users\wajordan\Desktop\2024_13_Annual_Review\PLOTS\AIRFOIL\ERR_OOA';
err_file = 'ERR_L2_CASE_3_IC-150.png';
ooa_file = 'OOA_L2_CASE_3_IC-150.png';


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