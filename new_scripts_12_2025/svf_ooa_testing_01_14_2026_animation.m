%% Parsing SVF data for SciTech (animations) (01/14/2026)
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

print_ERR=false;
print_OOA=false;
target_folder = 'C:\Users\wajordan\Desktop\';
err_file = '';
ooa_file = '';

%% for use with 'parse_and_plot_new2'
DATA_DIR='C:\Users\wajordan\Desktop\';
foldernames1 = {};
foldernames1 = [foldernames1,'SVF_2026-01-14_18.52.03_CURVED_P4_q6_q6_nohoex_TT_subsonic'];% 1
foldernames1 = [foldernames1,'SVF_2026-01-14_19.01.12_CURVED_P4_q6_q6_nohoex_FT_subsonic'];% 2

foldernames1 = [foldernames1,'SVF_2026-01-14_19.06.07_CURVED_P4_q6_q6_nohoex_FF_subsonic'];% 3
foldernames1 = [foldernames1,'SVF_2026-01-14_19.15.26_LINEAR_P4_q6_q6_nohoex_FF_subsonic'];% 4
foldernames1 = [foldernames1,'SVF_2026-01-14_19.21.02_LINEAR_P4_q6_q6_nohoex_TT_subsonic'];% 5

foldernames1 = [foldernames1,'SVF_2026-01-14_19.30.12_CURVED_P4_q6_q6_hoex_FF_subsonic'];                   % 6
foldernames1 = [foldernames1,'SVF_2026-01-14_19.38.09_CURVED_P4_q6_q6_hoex_FT_subsonic'];                   % 7

foldernames1 = [foldernames1,'SVF_2026-01-14_19.44.36_CURVED_P4_q6_q6_hoex_TT_subsonic'];                   % 8
foldernames1 = [foldernames1,'SVF_2026-01-14_21.31.31_CURVED_p4_q6_q6_hoex_TT_subsonic_geo4_10subiter'];    % 9

foldernames1 = [foldernames1,'SVF_2026-01-14_19.53.09_CURVED_p4_q6_q6_hoex_TT_BC_subsonic'];                % 10
foldernames1 = [foldernames1,'SVF_2026-01-15_01.12.50_CURVED_p4_q6_q6_hoex_TT_BC_subsonic_10subiter'];         % 11
foldernames1 = [foldernames1,'SVF_2026-01-14_20.42.31_CURVED_p4_q6_q6_hoex_TT_BC_subsonic_geo4'];           % 12
foldernames1 = [foldernames1,'SVF_2026-01-14_20.57.18_CURVED_p4_q6_q6_hoex_TT_BC_subsonic_geo4_10subiter']; % 13

foldernames1 = [foldernames1,'SVF_2026-01-14_20.14.59_CURVED_p4_q6_q6_hoex_TT_BC_supersonic'];      % 14
foldernames1 = [foldernames1,'SVF_2026-01-14_20.34.54_CURVED_p4_q6_q6_hoex_TT_BC_supersonic_geo4']; % 15
foldernames1 = [foldernames1,'SVF_2026-01-14_20.21.52_CURVED_p4_q6_q6_hoex_TT_supersonic'];         % 16
foldernames1 = [foldernames1,'SVF_2026-01-14_20.30.48_CURVED_p4_q6_q6_hoex_TT_supersonic_geo4'];    % 17
foldernames1 = cellfun(@(str_b)strcat(DATA_DIR,str_b),foldernames1,UniformOutput=false);

print_ERR=true;
print_OOA=true;
target_folder = 'C:\Users\wajordan\Desktop\';
% err_file = 'ERR_SUBSONIC_ALL_VARS_LINEAR_VS_CURVED.png';
% ooa_file = 'OOA_SUBSONIC_ALL_VARS_LINEAR_VS_CURVED.png';

err_file = 'ERR_SUBSONIC_ALL_VARS_EFFECT_OF_BCS.png';
ooa_file = 'OOA_SUBSONIC_ALL_VARS_EFFECT_OF_BCS.png';

% err_file = 'ERR_SUBSONIC_ALL_VARS_EFFECT_OF_B_ORDER.png';
% ooa_file = 'OOA_SUBSONIC_ALL_VARS_EFFECT_OF_B_ORDER.png';

% foldernames = foldernames1([9,9,13]);
foldernames = foldernames1([6,6,9]);
foldernames = foldernames1([11,11,13]);
var_select    = [ 3, 4, 4 ];
var_mask      = {[ 1, 1, 1, 1 ]};
norm_select   = [1];
layer_select  = {[]};
line_fmt      = { '-', ':','--' };
color_spec    = {lines(4)};
% tag_fmt       = { '', '(linear)', '(curved)' };
tag_fmt       = { '', '(curved)', '(curved w/ BC)' };
% tag_fmt       = { '', '(curved P2)', '(curved P4)' };
legend_flag   = true;



post_plot_commands = {"set(hfig1.Children(4),'Ylim',[1e-14,1e-2]);",...
                      "yticks(hfig1.Children(4),10.^(-14:2:-2));",  ...
                      "set(hfig1.Children(4),'Xlim',[8,256])",    ...
                      "set(hfig1.Children(3),'Location','southwest');",...
                      "set(hfig1.Children(2),'Ylim',[0,6]);",...
                      "set(hfig1.Children(2),'Xlim',[8,256])",...
                      "set(hfig1.Children(1),'Visible','off');"};

new = true;
for i = 150
    iter_select   = {[i],[i],[i]};
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

    label_iteration(hfig1.Children(4),i);
    label_iteration(hfig1.Children(2),i);
    cellfun(@eval,post_plot_commands)
    if (print_ERR)
        exportgraphics(hfig1.Children(4),fullfile(target_folder,err_file),'Resolution',600,Append=~new)
        % new = false;
    end
    if (print_OOA)
        exportgraphics(hfig1.Children(2),fullfile(target_folder,ooa_file),'Resolution',600,Append=~new)
        % new = false;
    end
end

function label_iteration(ax,iter)

xpos = 0.4;
ypos = 0.97;
txt = sprintf('Iterative Correction: %d',iter);

text(ax,xpos,ypos,txt,'Interpreter','latex','Units','Normalized',...
                      'VerticalAlignment','top')
end

function plot_orders(ax)

xpos = 0.4;
ypos = 0.97;
for i = 1:numel(ypos)
    plot([xlim],0.01*[xlim].^-5,'k--')
    txt = sprintf('%d',iter);
    text(ax,xpos,ypos,txt,'Interpreter','latex','Units','Normalized',...
                          'VerticalAlignment','top')
end
end