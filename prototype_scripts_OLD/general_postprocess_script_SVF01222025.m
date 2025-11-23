%% post-processing script
clc; clear; close all;
%%
dim   = 2;
r_fac = 2;
foldernames = {            'C:\Users\wajordan\Desktop\Desktop_TEMP\RESULTS\VORTEX\CASE_1\SVF_CASE1_2024-12-16_02.16.04'};
foldernames = [foldernames,'C:\Users\wajordan\Desktop\Desktop_TEMP\RESULTS\VORTEX\CASE_2\SVF_CASE2_2024-12-16_02.15.26'];
foldernames = [foldernames,'C:\Users\wajordan\Desktop\Desktop_TEMP\RESULTS\VORTEX\CASE_3\SVF_CASE3_2024-12-16_02.16.53'];
foldernames = [foldernames,'C:\Users\wajordan\Desktop\Desktop_TEMP\RESULTS\VORTEX\CASE_4a\SVF_CASE4a_2024-12-16_10.37.37'];
foldernames = [foldernames,'C:\Users\wajordan\Desktop\Desktop_TEMP\RESULTS\VORTEX\CASE_4b\SVF_CASE4b_2024-12-16_02.17.41'];

% foldernames = foldernames([5]);
% var_select    = [ 3, 4, 4, 4 ];
% var_mask      = {[ 1, 1, 1, 1 ]};
% norm_select   = [1];
% iter_select   = {[],[0],[1],[150]};
% layer_select  = {[]};
% line_fmt      = { '-', '--', '-.', ':' };
% color_spec    = {turbo(1),jet(1),hsv(1),lines(1)};
% % color_spec    = {lines(4)};
% tag_fmt       = { '' };
% legend_flag   = false;

foldernames = foldernames([2,4,5]);
var_select    = [ 3, 4, 4 ];
var_mask      = {[ 1, 1, 1, 1 ]};
norm_select   = [1];
iter_select   = {[],[50],[50]};
layer_select  = {[]};
line_fmt      = { '-', '--', '-.' };
color_spec    = {turbo(1),jet(1),hsv(1)};
% color_spec    = {lines(4)};
tag_fmt       = { '' };
legend_flag   = false;

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

