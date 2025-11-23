%% post-processing script
clc; clear; close all;
%%
dim   = 2;
r_fac = 2;

%% for use with 'parse_and_plot_new2'
% foldernames1 = {             'C:\Users\wajordan\Desktop\RESULTS\AIRFOIL\CASE_1\JOUKOWSKI_C_GRID_curved_bc_constraints_IC150_Iter150_alpha_0_UR0.5_2024-12-13_10.04.45'}; % 1
% foldernames1 = [foldernames1,'C:\Users\wajordan\Desktop\RESULTS\AIRFOIL\CASE_2\JOUKOWSKI_C_GRID_curved_bc_constraints_IC_full_alpha_1_UR0.5_2024-12-14_11.05.56'];       % 2
% foldernames1 = [foldernames1,'C:\Users\wajordan\Desktop\RESULTS\AIRFOIL\CASE_3\JOUKOWSKI_C_GRID_curved_2024-12-16_12.35.28'];
% foldernames1 = [foldernames1,'C:\Users\wajordan\Desktop\RESULTS\AIRFOIL\CASE_1\JOUKOWSKI_C_GRID_curved_bc_constraints_IC100_Iter100_alpha_0_UR0.5_2024-12-12_12.31.32'];



% foldernames1 = {             'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2025-02-12_09.01.33_OLD_STENCILS'};                 % 1
% foldernames1 = [foldernames1,'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2_2025-02-12_09.01.45_NEW_STENCILS'];               % 2
% 
% foldernames1 = [foldernames1,'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2025-02-13_10.12.07_OLD_STENCILS_rec=true'];        % 3
% foldernames1 = [foldernames1,'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2_2025-02-13_10.12.20_NEW_STENCILS_rec=true'];      % 4
% 
% foldernames1 = [foldernames1,'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2025-02-14_08.43.08_OLD_STENCILS_limiters_on'];     % 5
% foldernames1 = [foldernames1,'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2_2025-02-14_08.42.49_NEW_STENCILS_limiters_on'];   % 6
% 
% foldernames1 = [foldernames1,'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2025-02-12_23.42.00_OLD_STENCILS_no_constraints'];  % 7
% foldernames1 = [foldernames1,'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2_2025-02-12_23.41.00_NEW_STENCILS_no_constraints'];% 8

foldernames1 = {             'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2025-02-28_09.01.38'};                 % 1
foldernames1 = [foldernames1,'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2025-02-28_11.01.29_unconstrained'];               % 2

%% effect of iterative corrections alpha=0
% foldernames = foldernames1([1,1]);
% var_select    = [ 3, 4 ];
% var_mask      = {[ 1, 1, 1, 1 ]};
% norm_select   = [2];
% iter_select   = {[],[1]};
% layer_select  = {[]};
% line_fmt      = { '-', '--' };
% % color_spec    = {turbo(1),jet(1),hsv(1),lines(1)};
% color_spec    = {lines(4)};
% tag_fmt       = { '' };
% legend_flag   = false;
% post_plot_commands = {"set(hfig1.Children(4),'Ylim',[1e-12,1e0]);",...
%                       "yticks(hfig1.Children(4),10.^(-12:1:0));",  ...
%                       "set(hfig1.Children(4),'Xlim',[5,5000])",    ...
%                       "set(hfig1.Children(3),'Location','southwest');",...
%                       "set(hfig1.Children(2),'Ylim',[0,6]);",...
%                       "set(hfig1.Children(2),'Xlim',[5,5000])",...
%                       "set(hfig1.Children(1),'Visible','off');"};
% print_ERR=false;
% print_OOA=false;
% err_file = 'ERR.png';
% ooa_file = 'OOA.png';

%% effect of iterative corrections alpha=1
% foldernames = foldernames1([2,2,2]);
% var_select    = [ 3, 4, 4 ];
% var_mask      = {[ 1, 1, 1, 1 ]};
% norm_select   = [2];
% iter_select   = {[],[0],[100]};
% layer_select  = {[]};
% line_fmt      = { '-', '--', ':' };
% % color_spec    = {turbo(1),jet(1),hsv(1),lines(1)};
% color_spec    = {lines(4)};
% tag_fmt       = { '' };
% legend_flag   = false;
% post_plot_commands = {"set(hfig1.Children(4),'Ylim',[1e-12,1e0]);",...
%                       "yticks(hfig1.Children(4),10.^(-12:1:0));",  ...
%                       "set(hfig1.Children(4),'Xlim',[5,5000])",    ...
%                       "set(hfig1.Children(3),'Location','southwest');",...
%                       "set(hfig1.Children(2),'Ylim',[0,6]);",...
%                       "set(hfig1.Children(2),'Xlim',[5,5000])",...
%                       "set(hfig1.Children(1),'Visible','off');"};
% print_ERR=false;
% print_OOA=false;
% err_file = 'ERR.png';
% ooa_file = 'OOA.png';

%% effect of iterative corrections alpha=5
% foldernames = foldernames1([3,3,3]);
% var_select    = [ 3, 4, 4 ];
% var_mask      = {[ 1, 1, 1, 1 ]};
% norm_select   = [2];
% iter_select   = {[],[0],[150]};
% layer_select  = {[]};
% line_fmt      = { '-', '--', ':' };
% % color_spec    = {turbo(1),jet(1),hsv(1),lines(1)};
% color_spec    = {lines(4)};
% tag_fmt       = { '', '(0 iter.)', '(150 iter.)' };
% legend_flag   = true;
% post_plot_commands = {"set(hfig1.Children(4),'Ylim',[1e-10,1e-3]);",...
%                       "yticks(hfig1.Children(4),10.^(-10:1:-3));",  ...
%                       "set(hfig1.Children(4),'Xlim',[10,1000])",    ...
%                       "set(hfig1.Children(3),'Location','southwest');",...
%                       "set(hfig1.Children(2),'Ylim',[0,5]);",...
%                       "set(hfig1.Children(2),'Xlim',[10,1000])",...
%                       "set(hfig1.Children(1),'Visible','off');"};
% print_ERR=false;
% print_OOA=false;



foldernames = foldernames1([1]);
var_select    = [ 1 ];
var_mask      = {[ 0, 0, 0, 1 ]};
norm_select   = [1];
iter_select   = {[]};
layer_select  = {[]};
line_fmt      = { '-' };
% color_spec    = {turbo(1),jet(1),hsv(1),lines(1)};
color_spec    = {lines(4)};
tag_fmt       = { '' };
legend_flag   = true;
post_plot_commands = {"set(hfig1.Children(3),'Location','southwest');",...,...
                      "set(hfig1.Children(1),'Visible','off');"};


% foldernames = foldernames1([2,1,1]);
% var_select    = [ 3, 4, 4 ];
% var_mask      = {[ 1, 1, 1, 1 ]};
% norm_select   = [1];
% iter_select   = {[],[0],[2]};
% layer_select  = {[]};
% line_fmt      = { '-', '--', ':' };
% % color_spec    = {turbo(1),jet(1),hsv(1),lines(1)};
% color_spec    = {lines(4)};
% tag_fmt       = { '', '(1)', '(2)' };
% legend_flag   = true;
% post_plot_commands = {"set(hfig1.Children(4),'Ylim',[1e-10,1e-3]);",...
%                       "yticks(hfig1.Children(4),10.^(-10:1:-3));",  ...
%                       "set(hfig1.Children(4),'Xlim',[10,1000])",    ...
%                       "set(hfig1.Children(3),'Location','southwest');",...
%                       "set(hfig1.Children(2),'Ylim',[0,5]);",...
%                       "set(hfig1.Children(2),'Xlim',[10,1000])",...
%                       "set(hfig1.Children(1),'Visible','off');"};
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