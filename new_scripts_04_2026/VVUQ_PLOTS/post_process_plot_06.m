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

folders = [21,21,21,20];
prim = {'primal','primal','ic','ic'};
geom = {'linear','curved','curved','curved'};
% geom = {'linear','curved','linear','curved'};
legend_labels1 = {'$C_L$ primal','$C_L$ primal (reconstructed)','$C_L$ corr. (P2)','$C_L$ corr. (P4)'};
legend_labels2 = {'$C_D$ primal','$C_D$ primal (reconstructed)','$C_D$ corr. (P2)','$C_D$ corr. (P4)'};
x_label1 = '$N_{cells}^{1/2}$';
x_label2 = x_label1;
y_label1_1 ='$C_L$ Error';
y_label1_2 ='$C_D$ Error';
y_label2_1 ='$C_L$ OOA';
y_label2_2 ='$C_D$ OOA';
var1 = {'CL','CL','CL','CL'};
var2 = {'CD','CD','CD','CD'};
% ics = {[0],[0],[0],[0]};
% linspec ={'-','-','--',':'};
% colors  = lines(4);
% markers = {'o','s','>','<'};
linspec ={'-',':','--',':'};
colors  = validatecolor({'red','red','blue','green'},'multiple');
markers = {'o','o','d','d'};
marker_size = [9,9, 9,9];

post_plot_commands1 = {"set(hfig.Children(4),'Ylim',[1e-8,1e0]);",...
                      "yticks(hfig.Children(4),10.^(-8:1:0));",  ...
                      "set(hfig.Children(4),'Xlim',[10,1000])",    ...
                      "set(hfig.Children(3),'Location','southwest');",...
                      "set(hfig.Children(2),'Ylim',[0,7]);",...
                      "set(hfig.Children(2),'Xlim',[10,1000]);",...
                      "set(hfig.Children(1),'Location','southwest');"};




folders2 = {ALL_DATA.DATA([21,21,20]).folder};
tag_fmt = { '', '(P2)', '(P4)'};
var_select    = [ 3, 4, 4 ];
var_mask      = {[ 1, 1, 1, 1 ]};
norm_select   = [1];
layer_select  = {[]};
line_fmt      = { '-', '--', ':' };
color_spec    = {lines(4)};
legend_flag   = true;
post_plot_commands2 = {"set(hfig2.Children(4),'Ylim',[1e-11,1e-3]);",...
                      "yticks(hfig2.Children(4),10.^(-11:1:-3));",  ...
                      "set(hfig2.Children(4),'Xlim',[10,1000])",    ...
                      "set(hfig2.Children(3),'Location','southwest');",...
                      "set(hfig2.Children(2),'Ylim',[0,5]);",...
                      "set(hfig2.Children(2),'Xlim',[10,1000])",...
                      "set(hfig2.Children(1),'Visible','off');"};

target_folder = 'C:\Users\wajordan\Desktop\Meeting_Plots\ASME_VVUQ_2026_plots\ERR_CL';
name ='CL_CD_error_anim_ic20_P2vP4';

filename1 = fullfile(target_folder,[name,'.gif']);
filename2 = fullfile(target_folder,[name,'_start.png']);
filename3 = fullfile(target_folder,[name,'_end.png']);
framerate = 2;
end_iter = 20;
for i = 0:end_iter
    ics = {[0],[0],[i],[i]};
    [N,E,OOA] = gather_error_plot_data(ALL_DATA,folders,geom,prim,var1,ics);
    hfig1 = plot_OOA(N,E,OOA,linspec,colors,markers,marker_size,x_label1,y_label1_1,x_label2,y_label2_1,legend_labels1,post_plot_commands1);
    label_iteration(hfig1.Children(4),i);
    label_iteration(hfig1.Children(2),i);
    % set(hfig1.Children(4),'PlotBoxAspectRatio',[1,1.5,1])
    exportgraphics(hfig1.Children(4),'tmp2.tiff','Resolution',600)

    [N,E,OOA] = gather_error_plot_data(ALL_DATA,folders,geom,prim,var2,ics);
    hfig1 = plot_OOA(N,E,OOA,linspec,colors,markers,marker_size,x_label1,y_label1_2,x_label2,y_label2_2,legend_labels2,post_plot_commands1);
    label_iteration(hfig1.Children(4),i);
    label_iteration(hfig1.Children(2),i);
    % set(hfig1.Children(4),'PlotBoxAspectRatio',[1,1.5,1])
    exportgraphics(hfig1.Children(4),'tmp3.tiff','Resolution',600)

    iter_select = {[],[i],[i]}; 
    [hfig2,DE_test] = parse_and_plot_new2(2,2,folders2,          ...
                                                          var_select,   ...
                                                          var_mask,     ...
                                                          norm_select,  ...
                                                          iter_select,  ...
                                                          layer_select, ...
                                                          tag_fmt,      ...
                                                          line_fmt,     ...
                                                          color_spec,   ...
                                                          legend_flag );
    cellfun(@eval,post_plot_commands2);
    label_iteration(hfig2.Children(4),i);
    label_iteration(hfig2.Children(2),i);
    % set(hfig1.Children(4),'PlotBoxAspectRatio',[1,1.5,1])
    exportgraphics(hfig2.Children(4),'tmp1.tiff','Resolution',600)
    I = imtile(["tmp1.tiff","tmp2.tiff","tmp3.tiff"],BackgroundColor="w",GridSize=[1,3]);
    delete('tmp1.tiff','tmp2.tiff','tmp3.tiff');

    [A,map] = rgb2ind(I,256);
    if i ==0
        imwrite(A,map,filename1,"gif",LoopCount=Inf,DelayTime=1/framerate)
        imwrite(A,map,filename2,"png");
    elseif i==end_iter
        imwrite(A,map,filename1,"gif",WriteMode="append",DelayTime=1/framerate)
        imwrite(A,map,filename3,"png");
    else
        imwrite(A,map,filename1,"gif",WriteMode="append",DelayTime=1/framerate)
    end
    hold off
    % exportgraphics(hfig.Children(4),filename,'Resolution',600,'Append',true)
    % drawnow
end


function label_iteration(ax,iter)
xpos = 0.4;
ypos = 0.97;
txt = sprintf('Iterative Correction: %d',iter);
text(ax,xpos,ypos,txt,'Interpreter','latex','Units','Normalized',...
                      'VerticalAlignment','top')
end