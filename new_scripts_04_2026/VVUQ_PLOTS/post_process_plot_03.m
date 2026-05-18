%% Parsing KT-airfoil data post-process: cweno reconstruction compare (05/15/2026)
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
DATA_DIR='C:\Users\wajordan\Desktop\git_MATLAB\SENSEI_POSTPROCESSING\new_scripts_04_2026\data_preprocess';
file_name = 'F_ALL_05_16_2026.mat';
load(fullfile(DATA_DIR,file_name));

% DATA_DIR='C:\Users\wajordan\Desktop\CASES\';
% foldernames1 = {};
% foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2026-01-15_02.00.10_CURVED_P4_q6_q6_hoex_TT_BC_VARLAY2_ITER2000_IC10'];% 1
% foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2026-01-15_03.10.03_CURVED_P4_q6_q6_hoex_TT_BC_VARLAY8_ITER2000_IC10'];% 2
% foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2026-02-13_12.01.06_ORDER_4_bc_IC_10_100_iter_GEO4'];% 3
% foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_limited_2026-02-13_12.00.09_ORDER_4_bc_IC_10_100_iter_GEO4'];% 4
% foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_04_01_2026-04-02_12.20.09'];        % 5
% foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_04_01_regress_2026-04-02_12.21.04'];% 6
% 
% foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_clustered_04_02_regress_2026-04-02_14.54.19'];% 7
% foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_clustered_04_02_2026-04-02_14.54.15'];% 8
% 
% foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_clustered_04_02_regress_2026-04-02_19.10.44'];% 9
% foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_clustered_04_02_2026-04-02_19.10.39'];% 10
% 
% foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_clustered_04_02_regress_2026-04-03_11.01.28'];% 11
% foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_clustered_04_02_2026-04-03_10.47.17'];% 12
% 
% 
% foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_v2_04_08_regress_2026-04-09_11.38.29']; % 13
% foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_v2_04_08_2026-04-09_11.38.25'        ]; % 14
% 
% foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_v3_04_09_regress_2026-04-09_18.12.43'        ]; % 15
% foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_v3_04_09_2026-04-09_18.12.39'        ]; % 16
% 
% foldernames1 = [foldernames1,  'JOUKOWSKI_C_GRID_curved_v3_04_09_regress_2026-04-10_11.40.43'      ]; % 17
% foldernames1 = [foldernames1,  'JOUKOWSKI_C_GRID_curved_v3_04_09_2026-04-10_11.40.38'      ]; % 18
% 
% foldernames1 = [foldernames1,  'REGRESS_4_JOUKOWSKI_C_GRID_curved_v3_04_09_regress_2026-05-13_12.36.01'      ]; % 19
% 
% foldernames1 = [foldernames1,  'ITER_0_JOUKOWSKI_C_GRID_curved_v3_04_09_regress_2026-05-14_20.02.21'      ]; % 20
% foldernames1 = [foldernames1,  'P2_ITER_0_JOUKOWSKI_C_GRID_curved_v3_04_09_regress_2026-05-15_11.04.23'      ]; % 21

folders = [17,17,17,17];
prim = {'primal','primal','ic','ic'};
geom = {'linear','curved','linear','curved'};
% legend_labels = {'$C_L$ Primal','$C_L$ Primal (reconstructed)','$C_L$ ETE (k-exact, reconstructed)','$C_L$ ETE (CWENO, reconstructed)'};
legend_labels = {'$C_L$ Primal','$C_L$ Primal (reconstructed)','$C_L$ ETE (k-exact)','$C_L$ ETE (k-exact, reconstructed)'};
x_label1 = '$N_{cells}^{1/2}$';
x_label2 = x_label1;
y_label1 ='$C_L$ Error';
y_label2 ='$C_L$ OOA';
var = {'CL','CL','CL','CL'};
% ics = {[0],[0],[0],[0]};
linspec ={'-','-',...
          '--','--'};
colors  = lines(4);
markers = {'^','v','o','s'};
marker_size = [9,9, 9,9];

post_plot_commands = {"set(hfig.Children(4),'Ylim',[1e-7,1e0]);",...
                      "yticks(hfig.Children(4),10.^(-7:1:0));",  ...
                      "set(hfig.Children(4),'Xlim',[10,1000])",    ...
                      "set(hfig.Children(3),'Location','southwest');",...
                      "set(hfig.Children(2),'Ylim',[0,7]);",...
                      "set(hfig.Children(2),'Xlim',[10,1000]);",...
                      "set(hfig.Children(1),'Location','southwest');"};


target_folder = 'C:\Users\wajordan\Desktop\Meeting_Plots\ASME_VVUQ_2026_plots\ERR_CL';
name ='CL_error_anim_ic20_linear_v_curved.mp4';
v = VideoWriter(fullfile(target_folder,name),"MPEG-4");
v.Quality = 100;
v.FrameRate = 1;
open(v);
for i = 0:20
    ics = {[0],[0],[i],[i]};
    [N,E,OOA] = gather_error_plot_data(ALL_DATA,folders,geom,prim,var,ics);
    hfig = plot_OOA(N,E,OOA,linspec,colors,markers,marker_size,x_label1,y_label1,x_label2,y_label2,legend_labels,post_plot_commands);
    label_iteration(hfig.Children(4),i);
    label_iteration(hfig.Children(2),i);
    cellfun(@eval,post_plot_commands)
    F = print('-RGBImage','-r300');
    % writeVideo(v,F)
    

    % extra frame in case matlab fucks it up
    if (i==0)
        writeVideo(v,F)
    end
    hold off
    drawnow
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