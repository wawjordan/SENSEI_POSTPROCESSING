%% post-processing script
clc; clear; close all;
%%
dim   = 2;
r_fac = 2;

%% for use with 'parse_and_plot_new2'
foldernames1 = {            'C:\Users\wajordan\Desktop\RESULTS\VORTEX\CASE_1\SVF_CASE1_2024-12-16_02.16.04'};   % 1
foldernames1 = [foldernames1,'C:\Users\wajordan\Desktop\RESULTS\VORTEX\CASE_2\SVF_CASE2_2024-12-16_02.15.26'];  % 2
foldernames1 = [foldernames1,'C:\Users\wajordan\Desktop\RESULTS\VORTEX\CASE_3\SVF_CASE3_2024-12-16_02.16.53'];  % 3
foldernames1 = [foldernames1,'C:\Users\wajordan\Desktop\RESULTS\VORTEX\CASE_4a\SVF_CASE4a_2024-12-16_10.37.37'];% 4
foldernames1 = [foldernames1,'C:\Users\wajordan\Desktop\RESULTS\VORTEX\CASE_4b\SVF_CASE4b_2024-12-16_02.17.41'];% 5


%% CASE 1 & CASE 2 (no iterative correction)
% foldernames = foldernames1([1,1,2]);
% var_select    = [ 3, 4, 4 ];
% var_mask      = {[ 1, 1, 1, 1 ]};
% norm_select   = [2];
% iter_select   = {[],[0],[0]};
% layer_select  = {[]};
% line_fmt      = { '-', '--', ':' };
% color_spec    = {lines(4)};
% tag_fmt       = { '','(linear)','(curved)' };
% legend_flag   = true;
% post_plot_commands = {"set(hfig1.Children(4),'Ylim',[1e-8,1e-2]);",...
%                       "yticks(hfig1.Children(4),10.^(-8:1:-2));",  ...
%                       "set(hfig1.Children(4),'Xlim',[5,1000])",    ...
%                       "set(hfig1.Children(3),'Location','southwest');",...
%                       "set(hfig1.Children(2),'Ylim',[1.8,2.2]);",...
%                       "set(hfig1.Children(2),'Xlim',[5,1000])",...
%                       "set(hfig1.Children(1),'Visible','off');"};
% print_ERR=true;
% print_OOA=true;
% target_folder = 'C:\Users\wajordan\Desktop\2024_13_Annual_Review\PLOTS\VORTEX\ERR_OOA';
% err_file = 'ERR_L2_CASE_1-2_No-IC.png';
% ooa_file = 'OOA_L2_CASE_1-2_No-IC.png';

%% CASE 1 & CASE 2 (10 iterative corrections)
% foldernames = foldernames1([1,1,2]);
% var_select    = [ 3, 4, 4 ];
% var_mask      = {[ 1, 1, 1, 1 ]};
% norm_select   = [2];
% iter_select   = {[],[10],[10]};
% layer_select  = {[]};
% line_fmt      = { '-', '--', ':' };
% color_spec    = {lines(4)};
% tag_fmt       = { '','(linear)','(curved)' };
% legend_flag   = true;
% post_plot_commands = {"set(hfig1.Children(4),'Ylim',[1e-10,1e-2]);",...
%                       "yticks(hfig1.Children(4),10.^(-10:1:-2));",  ...
%                       "set(hfig1.Children(4),'Xlim',[5,1000])",    ...
%                       "set(hfig1.Children(3),'Location','southwest');",...
%                       "set(hfig1.Children(2),'Ylim',[1,3]);",...
%                       "set(hfig1.Children(2),'Xlim',[5,1000])",...
%                       "set(hfig1.Children(1),'Visible','off');"};
% print_ERR=true;
% print_OOA=true;
% target_folder = 'C:\Users\wajordan\Desktop\2024_13_Annual_Review\PLOTS\VORTEX\ERR_OOA';
% err_file = 'ERR_L2_CASE_1-2_IC-10.png';
% ooa_file = 'OOA_L2_CASE_1-2_IC-10.png';


%% CASE 1 & CASE 2 (iterative correction to convergence)
% foldernames = foldernames1([1,1,2]);
% var_select    = [ 3, 4, 4 ];
% var_mask      = {[ 1, 1, 1, 1 ]};
% norm_select   = [2];
% iter_select   = {[],[],[]};
% layer_select  = {[]};
% line_fmt      = { '-', '--', ':' };
% color_spec    = {lines(4)};
% tag_fmt       = { '','(linear)','(curved)' };
% legend_flag   = true;
% post_plot_commands = {"set(hfig1.Children(4),'Ylim',[1e-14,1e-2]);",...
%                       "yticks(hfig1.Children(4),10.^(-14:2:-2));",  ...
%                       "set(hfig1.Children(4),'Xlim',[5,1000])",    ...
%                       "set(hfig1.Children(3),'Location','southwest');",...
%                       "set(hfig1.Children(2),'Ylim',[0,6]);",...
%                       "set(hfig1.Children(2),'Xlim',[5,1000])",...
%                       "set(hfig1.Children(1),'Visible','off');"};
% print_ERR=true;
% print_OOA=true;
% target_folder = 'C:\Users\wajordan\Desktop\2024_13_Annual_Review\PLOTS\VORTEX\ERR_OOA';
% err_file = 'ERR_L2_CASE_1-2_IC-to-convergence.png';
% ooa_file = 'OOA_L2_CASE_1-2_IC-to-convergence.png';

%% CASE 2 vs 4a (0 Iterative Corrections)
% foldernames = foldernames1([2,2,4]);
% var_select    = [ 3, 4, 4 ];
% var_mask      = {[ 1, 1, 1, 1 ]};
% norm_select   = [2];
% iter_select   = {[],[0],[0]};
% layer_select  = {[]};
% line_fmt      = { '-', '--', ':' };
% color_spec    = {lines(4)};
% tag_fmt       = { '','(unconstrained)','(constrained)' };
% legend_flag   = true;
% post_plot_commands = {"set(hfig1.Children(4),'Ylim',[1e-9,1e-2]);",...
%                       "yticks(hfig1.Children(4),10.^(-9:1:-2));",  ...
%                       "set(hfig1.Children(4),'Xlim',[5,1000])",    ...
%                       "set(hfig1.Children(3),'Location','southwest');",...
%                       "set(hfig1.Children(2),'Ylim',[1.8,2.2]);",...
%                       "set(hfig1.Children(2),'Xlim',[5,1000])",...
%                       "set(hfig1.Children(1),'Visible','off');"};
% print_ERR=true;
% print_OOA=true;
% target_folder = 'C:\Users\wajordan\Desktop\2024_13_Annual_Review\PLOTS\VORTEX\ERR_OOA';
% err_file = 'ERR_L2_CASE_2-4_No-IC.png';
% ooa_file = 'OOA_L2_CASE_2-4_No-IC.png';

%% CASE 2 vs 4a (10 Iterative Corrections)
% foldernames = foldernames1([2,2,4]);
% var_select    = [ 3, 4, 4 ];
% var_mask      = {[ 1, 1, 1, 1 ]};
% norm_select   = [2];
% iter_select   = {[],[10],[10]};
% layer_select  = {[]};
% line_fmt      = { '-', '--', ':' };
% color_spec    = {lines(4)};
% tag_fmt       = { '','(unconstrained)','(constrained)' };
% legend_flag   = true;
% post_plot_commands = {"set(hfig1.Children(4),'Ylim',[1e-11,1e-2]);",...
%                       "yticks(hfig1.Children(4),10.^(-11:1:-2));",  ...
%                       "set(hfig1.Children(4),'Xlim',[5,1000])",    ...
%                       "set(hfig1.Children(3),'Location','southwest');",...
%                       "set(hfig1.Children(2),'Ylim',[0,6]);",...
%                       "set(hfig1.Children(2),'Xlim',[5,1000])",...
%                       "set(hfig1.Children(1),'Visible','off');"};
% print_ERR=true;
% print_OOA=true;
% target_folder = 'C:\Users\wajordan\Desktop\2024_13_Annual_Review\PLOTS\VORTEX\ERR_OOA';
% err_file = 'ERR_L2_CASE_2-4_IC-10.png';
% ooa_file = 'OOA_L2_CASE_2-4_IC-10.png';


%% CASE 2 vs 4a (20 Iterative Corrections)
foldernames = foldernames1([2,2,4]);
var_select    = [ 3, 4, 4 ];
var_mask      = {[ 1, 1, 1, 1 ]};
norm_select   = [2];
iter_select   = {[],[50],[50]};
layer_select  = {[]};
line_fmt      = { '-', '--', ':' };
color_spec    = {lines(4)};
tag_fmt       = { '','(unconstrained)','(constrained)' };
legend_flag   = true;
post_plot_commands = {"set(hfig1.Children(4),'Ylim',[1e-16,1e-2]);",...
                      "yticks(hfig1.Children(4),10.^(-16:2:-2));",  ...
                      "set(hfig1.Children(4),'Xlim',[5,1000])",    ...
                      "set(hfig1.Children(3),'Location','southwest');",...
                      "set(hfig1.Children(2),'Ylim',[0,6]);",...
                      "set(hfig1.Children(2),'Xlim',[5,1000])",...
                      "set(hfig1.Children(1),'Visible','off');"};
print_ERR=true;
print_OOA=true;
target_folder = 'C:\Users\wajordan\Desktop\2024_13_Annual_Review\PLOTS\VORTEX\ERR_OOA';
err_file = 'ERR_SEQUENCE_L2_CASE_2-4_IC-50.png';
ooa_file = 'OOA_SEQUENCE_L2_CASE_2-4_IC-50.png';


%%
[hfig1,DE_test] = parse_and_plot_new2( dim, r_fac, foldernames,         ...
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
% set(hfig1.Children(4),'Xlim',[1,1000])
% set(hfig1.Children(3),'Visible','off');
% set(hfig1.Children(3),'Location','southwest');

% set(hfig1.Children(2),'Ylim',[0,6])
% yticks(hfig1.Children(2),0:6)
% set(hfig1.Children(2),'Xlim',[1,1000])
% set(hfig1.Children(1),'Visible','off');
% set(hfig1.Children(1),'Location','southwest');


% set(hfig1.Children(3),'XlimMode','auto')
% set(hfig1.Children(2),'XLimMode','auto')
% hfig1.Children(2).Legend.Location = 'southwest';



% target_folder = fullfile('E:\wajordan\Desktop_Overflow/',string(datetime(datetime(),'Format','yyyy-MM-dd_HH.mm.ss')));