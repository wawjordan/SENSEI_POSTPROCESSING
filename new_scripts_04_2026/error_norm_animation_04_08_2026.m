%% Parsing SVF data for SciTech (animations) (04/08/2026)
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
foldernames1 = {             'JOUKOWSKI_C_GRID_curved_04_01_2026-04-02_12.20.09'                  }; % 1
foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_04_01_regress_2026-04-02_12.21.04'          ]; % 2
foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_clustered_04_02_regress_2026-04-02_14.54.19']; % 3
foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_clustered_04_02_2026-04-02_14.54.15'        ]; % 4
foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_clustered_04_02_regress_2026-04-02_19.10.44']; % 5
foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_clustered_04_02_2026-04-02_19.10.39'        ]; % 6
foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_clustered_04_02_regress_2026-04-03_11.01.28']; % 7
foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_clustered_04_02_2026-04-03_10.47.17'        ]; % 8
foldernames1 = cellfun(@(str_b)strcat(DATA_DIR,str_b),foldernames1,UniformOutput=false);

print_ERR=true;
print_OOA=true;

% foldernames = foldernames1([8,7,8]);
foldernames = foldernames1([6,5,6]);
var_select    = [ 3, 4, 4 ];
var_mask      = {[ 1, 1, 1, 1 ]};
norm_select   = [1];
layer_select  = {[]};
line_fmt      = { '-', ':','--' };
color_spec    = {lines(4)};
tag_fmt       = { '', '(K-exact)', '(CWENO)' };
legend_flag   = true;



post_plot_commands = {"set(hfig1.Children(4),'Ylim',[1e-14,1e-2]);",...
                      "yticks(hfig1.Children(4),10.^(-14:2:-2));",  ...
                      "set(hfig1.Children(4),'Xlim',[10,1000])",    ...
                      "set(hfig1.Children(3),'Location','southwest');",...
                      "set(hfig1.Children(2),'Ylim',[0,6]);",...
                      "set(hfig1.Children(2),'Xlim',[10,1000])",...
                      "set(hfig1.Children(1),'Visible','off');"};
target_folder = 'C:\Users\wajordan\Desktop\Meeting_Plots\research_meeting_04_08_2026';
v = VideoWriter(fullfile(target_folder,'L1_error_anim.mp4'),"MPEG-4");
v.Quality = 100;
v.FrameRate = 10;
open(v);
for i = 0:200
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
    F = print('-RGBImage','-r300');
    writeVideo(v,F)

    % extra frame in case matlab fucks it up
    if (i==0)
        writeVideo(v,F)
    end
end
% extra frame in case matlab fucks it up
writeVideo(v,F)
close(v);

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