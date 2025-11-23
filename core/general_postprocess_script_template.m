%% post-processing script template
clc; clear; close all;
%%
% dimension - used for determining cell count for plotting error norms
dim   = 2;

% refinement factor between successive grids for the refinement study
r_fac = 2;

%% for use with 'parse_and_plot_new2'
% foldernames = {'<output-folder1-from-copy_results_to_desktop.sh>'};
% foldernames = [foldernames,'<output-folder1-from-copy_results_to_desktop.sh>'];


%% example inputs

% select which folders to compare (can just be 1)
foldernames = foldernames([1]);

% select the variables to read
%   1 - reconstruction error (primitive variables)
%   2 - error in TE estimate
%   3 - primal solution DE (primitive variables)
%   4 - corrected (ETE) solution DE (primitive variables)
%   5 - primal solution DE (conserved variables)
%   6 - corrected (ETE) solution DE (conserved variables)
%   each entry corresponds to a selected folder
var_select    = [ 3 ];

% select which variables to plot (2D currently)
%   [rho,u,v,P] if primitive variables
%   [ rho, rho u, rho v, rho e_t ] if conserved variables
%   set mask entries to 1 to include in plot, 0 to exclude
%   each mask [] in the cell array {} corresponds to a selected folder
var_mask      = {[ 1, 0, 0, 0 ]};

% select which norms to plot
%  1 - Discrete L_1 norm
%  2 - Discrete L_2 norm
%  3 - Discrete L_infinity norm
%   each entry corresponds to a selected folder
norm_select   = [ 3 ];

% select which iterative corrections to print
% (set as empty brackets if no iterative corrections)
%   each range [] in the cell array {} corresponds to a selected folder
iter_select   = {[]};

% if excluding layers from error norm calculations, how many?
% (set as empty brackets if not using)
%   each range [] in the cell array {} corresponds to a selected folder
layer_select  = {[]};

% plotting format for the lines
%   each entry in the cell array {} corresponds to a selected folder
line_fmt      = { '-' };

% plotting color for the lines (this dimension will correspond to the
%                               different variables)
%   each entry in the cell array {} corresponds to a selected folder
color_spec    = {lines(4)};

% additional text for displaying in the legend
%   each entry in the cell array {} corresponds to a selected folder
tag_fmt       = { '' };

% flag to denote if the legend should be printed for each entry (false),
% or if a summarized legend should be displayed (true)
legend_flag   = false;

% post plotting commands (details will vary based on the output)
post_plot_commands = {"set(hfig1.Children(4),'Ylim',[1e-10,1e-3]);",...
                      "yticks(hfig1.Children(4),10.^(-10:1:-3));",  ...
                      "set(hfig1.Children(4),'Xlim',[10,1000])",    ...
                      "set(hfig1.Children(3),'Location','southwest');",...
                      "set(hfig1.Children(2),'Ylim',[0,5]);",...
                      "set(hfig1.Children(2),'Xlim',[10,1000])",...
                      "set(hfig1.Children(1),'Visible','off');"};


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
cellfun(@eval,post_plot_commands)