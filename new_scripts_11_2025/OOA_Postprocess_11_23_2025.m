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
DATA_DIR='';
foldernames1 = {             ''};% 1
foldernames1 = [foldernames1,''];% 2

foldernames = foldernames1([2,1,1]);
var_select    = [ 3, 4, 4 ];
var_mask      = {[ 1, 1, 1, 1 ]};
norm_select   = [1];
iter_select   = {[],[0],[2]};
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
target_folder = '';
err_file = 'ERR.png';
ooa_file = 'OOA.png';

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