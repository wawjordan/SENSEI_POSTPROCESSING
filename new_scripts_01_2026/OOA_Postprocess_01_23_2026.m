%% Parsing KT-airfoil data (01/23/2026)
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
DATA_DIR='C:\Users\Will\Desktop\MY_CASES\';
foldernames1 = {             'cts_new_2026-01-23_15.56.26'};% 1
foldernames1 = [foldernames1,'cts_new_2026-01-23_16.38.08'];% 2

foldernames1 = cellfun(@(str_b)strcat(DATA_DIR,str_b),foldernames1,UniformOutput=false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
foldernames = foldernames1([1,1,2]);
var_select    = [ 3, 4, 4 ];
var_mask      = {[ 1, 1, 1, 0, 1, 1 ]};
norm_select   = [1];
iter_select   = {[]};
layer_select  = {[]};
line_fmt      = { '-', '--', ':' };
color_spec    = {hsv(7)};
tag_fmt       = { '' };
legend_flag   = true;
post_plot_commands = {"set(hfig1.Children(4),'Ylim',[1e-11,1e-3]);",...
                      "yticks(hfig1.Children(4),10.^(-11:1:-3));",  ...
                      "set(hfig1.Children(4),'Xlim',[10,1000])",    ...
                      "set(hfig1.Children(3),'Location','southwest');",...
                      "set(hfig1.Children(2),'Ylim',[0,5]);",...
                      "set(hfig1.Children(2),'Xlim',[10,1000])",...
                      "set(hfig1.Children(1),'Visible','off');"};
print_ERR=false;
print_OOA=false;
target_folder = 'C:\Users\Will\Desktop';
err_file = 'ERR.png';
ooa_file = 'OOA.png';

[hfig1,DE_test] = parse_and_plot_new3(dim,r_fac, foldernames,          ...
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