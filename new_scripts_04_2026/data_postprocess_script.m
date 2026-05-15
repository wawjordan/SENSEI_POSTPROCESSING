%% Parsing KT-airfoil data post-process (05/14/2026)
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
file_name = 'F_17_18.mat';
load(fullfile(DATA_DIR,file_name));



% folders = [1,1,1,1];
% prim = {'primal','primal','ic','ic'};
% geom = {'linear','curved','linear','curved'};
% legend_labels = {'Primal','Primal (reconstructed)', 'ETE','ETE (reconstructed)'};
% markers = {'m-o','r-o','b:s','g--d'};
% var = {'CL','CL','CL','CL'};
% ics = {[0],[0],[200],[200]};


% folders = [1,1,1,1,1,1,1,1];
% prim = {'primal','primal','primal','primal','ic','ic','ic','ic'};
% geom = {'linear','curved','linear','curved','linear','curved','linear','curved'};
% legend_labels = {'$C_L$ Primal','$C_L$ Primal (reconstructed)','$C_D$ Primal','$C_D$ Primal (reconstructed)',...
%                  '$C_L$ ETE','$C_L$ ETE (reconstructed)','$C_D$ ETE','$C_D$ ETE (reconstructed)'};
% markers = {'r-^','r--v','b-s','b--s','m-^','m--v','g-s','g--s'};
% var = {'CL','CL','CD','CD','CL','CL','CD','CD'};
% ics = {[0],[0],[0],[0],[0],[0],[0],[0]};


% folders = [1,1,1,1,...
%            1,1,1,1,...
%            2,2,2,2];
% prim = {'primal','primal','primal','primal',...
%         'ic','ic','ic','ic',...
%         'ic','ic','ic','ic'};
% geom = {'linear','curved','linear','curved',...
%         'linear','curved','linear','curved',...
%         'linear','curved','linear','curved'};
% legend_labels = {'$C_L$ Primal','$C_L$ Primal (reconstructed)','$C_D$ Primal','$C_D$ Primal (reconstructed)',...
%                  '$C_L$ ETE (k-exact)','$C_L$ ETE (k-exact, reconstructed)','$C_D$ ETE (k-exact)','$C_D$ ETE (k-exact, reconstructed)',...
%                  '$C_L$ ETE (CWENO)','$C_L$ ETE (CWENO, reconstructed)','$C_D$ ETE (CWENO)','$C_D$ ETE (CWENO, reconstructed)'};
% markers = {'r-^','r--^','b-s','b--s',...
%            'm-^','m--^','g-p','g--p',...
%            'm-o','m--o','g-p','g--p'};
% var = {'CL','CL','CD','CD',...
%        'CL','CL','CD','CD',...
%        'CL','CL','CD','CD'};
% ics = {[0],[0],[0],[0],...
%        [0],[0],[0],[0],...
%        [0],[0],[0],[0]};

folders = [1,1,...
           1,1,...
           2,2];
prim = {'primal','primal',...
        'ic','ic',...
        'ic','ic'};
geom = {'linear','curved',...
        'linear','curved',...
        'linear','curved'};
legend_labels = {'$C_L$ Primal','$C_L$ Primal (reconstructed)',...
                 '$C_L$ ETE (k-exact)','$C_L$ ETE (k-exact, reconstructed)',...
                 '$C_L$ ETE (CWENO)','$C_L$ ETE (CWENO, reconstructed)'};
linspec ={'-','--',...
          '-','--',...
          '-','--'};
colors  = reshape([lines(3),lines(3)].',3,[]).';
markers = {'^','^',...
           '^','^',...
           'o','o'};
marker_size = [9,9, 9,9, 9,9];
var = {'CL','CL',...
       'CL','CL',...
       'CL','CL'};
ics = {[0],[0],...
       [0],[0],...
       [0],[0]};


N = retrieve_variable( ALL_DATA, folders(1), geom{1}, '', 'N' );
E = {};
OOA = {};
for i = 1:numel(folders)
    tmp = abs(retrieve_variable_error( ALL_DATA, folders(i), geom{i}, prim{i}, var{i} ));
    tmp = tmp(ics{i}(:)+1,:);
    E = [E,tmp];
    tmp = retrieve_variable_ooa( ALL_DATA, folders(i), geom{i}, prim{i}, var{i} );
    tmp = tmp(ics{i}(:)+1,:);
    OOA = [OOA,tmp];
end
post_plot_commands = {"set(hfig.Children(4),'Ylim',[1e-6,1e0]);",...
                      "yticks(hfig.Children(4),10.^(-6:1:0));",  ...
                      "set(hfig.Children(4),'Xlim',[10,1000])",    ...
                      "set(hfig.Children(3),'Location','southwest');",...
                      "set(hfig.Children(2),'Ylim',[0,5]);",...
                      "set(hfig.Children(2),'Xlim',[10,1000]);",...
                      "set(hfig.Children(1),'Location','southwest');"};
x_label1 = '$N_{cells}^{1/2}$';
x_label2 = x_label1;
y_label1 ='Error';
y_label2 ='OOA';
% hfig = plot_OOA(N,E,OOA,markers,x_label1,y_label1,x_label2,y_label2,legend_labels,post_plot_commands);
hfig = plot_OOA(N,E,OOA,linspec,colors,markers,marker_size,x_label1,y_label1,x_label2,y_label2,legend_labels,post_plot_commands);
hold off

% ind = 4;
% airfoil = make_airfoil(ALL_DATA);
% theta_interval = retrieve_variable( ALL_DATA, 1, geom, '', 'theta_interval' ); theta_interval = theta_interval{ind};
% x = retrieve_variable( ALL_DATA, 1, geom, '', 'x' ); x = x{ind};
% y = retrieve_variable( ALL_DATA, 1, geom, '', 'y' ); y = y{ind};
% exact_sim_cp1   = retrieve_variable( ALL_DATA, 1, 'linear', 'exact_sim', 'CP' );
% exact_sim_cp1 = exact_sim_cp1{ind}(:);
% exact_sim_cp2 = [exact_sim_cp1.';exact_sim_cp1.']; exact_sim_cp2 = exact_sim_cp2(:);
% 
% plot(theta_interval,exact_sim_cp2,'b');
% hold on
% airfoil.plot_piecewise_constant_data_theta(x,y,exact_sim_cp1,true,'r:');
% fplot(@(theta)airfoil.CP(theta),[0,2*pi],'r')
% hold off;

function variable = retrieve_variable_error( ALL_DATA, folder_num, geom, prim, var )
    variable = retrieve_variable( ALL_DATA, folder_num, geom, prim, var );
    airfoil = make_airfoil(ALL_DATA);
    variable = variable-airfoil.(var);
end
function variable = retrieve_variable_ooa( ALL_DATA, folder_num, geom, prim, var )
    error = abs(retrieve_variable_error( ALL_DATA, folder_num, geom, prim, var ));
    N = retrieve_variable( ALL_DATA, 1, geom, '', 'N' );
    variable = nan*error;
    r_fac = 2;
    for i = 2:numel(N)
        variable(:,i) = log(error(:,i-1) ./ error(:,i))./log(r_fac);
    end
end
function variable = retrieve_variable( ALL_DATA, folder_num, geom, prim, var )
inputs = ALL_DATA.inputs;
if strcmp(var,'N')
    variable = [ALL_DATA.DATA(folder_num).(geom).F(:).N];
    variable = variable/inputs.nskip;
elseif strcmp(var,'XC')
    variable = {ALL_DATA.DATA(folder_num).(geom).F(:).XC};
elseif strcmp(var,'xc')||strcmp(var,'yc')||...
       strcmp(var,'x') ||strcmp(var,'y')||...
       strcmp(var,'x0')||strcmp(var,'y0')
    N = numel(ALL_DATA.DATA(folder_num).(geom).G);
    variable = cell(N,1);
    for i = 1:N
        variable{i} = ALL_DATA.DATA(folder_num).(geom).G(i).(var);
    end
elseif strcmpi(var,'theta_interval')
    airfoil = make_airfoil(ALL_DATA);
    N = numel(ALL_DATA.DATA(folder_num).(geom).G);
    variable = cell(N,1);
    for i = 1:N
        X = ALL_DATA.DATA(folder_num).(geom).G(i).x;
        Y = ALL_DATA.DATA(folder_num).(geom).G(i).y;
        variable{i} = airfoil.get_theta_from_coords_piecewise_constant(X,Y,true);
    end
elseif strcmpi(var,'theta_center')
    airfoil = make_airfoil(ALL_DATA);
    N = numel(ALL_DATA.DATA(folder_num).(geom).G);
    variable = cell(N,1);
    for i = 1:N
        X = ALL_DATA.DATA(folder_num).(geom).G(i).xc;
        Y = ALL_DATA.DATA(folder_num).(geom).G(i).yc;
        tmp = airfoil.get_theta_from_coords(X,Y,true);
        variable{i} = tmp(1:2:end);
    end
elseif strcmp(var,'CP')||strcmp(var,'P')
    N = numel(ALL_DATA.DATA(folder_num).(geom).F);
    variable = cell(N,1);
    for i = 1:N
        variable{i} = ALL_DATA.DATA(folder_num).(geom).F(i).([prim,'_',var]);
    end
else
    variable = [ALL_DATA.DATA(folder_num).(geom).H(:).([prim,'_',var])];
end
end

function airfoil = make_airfoil(ALL_DATA)
inputs = ALL_DATA.inputs;
airfoil        = kt_airfoil( inputs.epsilon, inputs.kappa, inputs.tau );
airfoil.vinf   = inputs.vinf;
airfoil.rhoinf = inputs.rhoinf;
airfoil.pinf   = inputs.pinf;
airfoil        = airfoil.set_alpha(inputs.alpha);
end

function [hfig,p] = plot_OOA(N,E,OOA,linspec,colors,markers,marker_size,x_label1,y_label1,x_label2,y_label2,legend_labels,post_plot_commands)
lim1 = 10^( floor( log10( N(1)   ) ) );
lim2 = 10^( ceil(  log10( N(end) ) ) );

hfig = stdplot(1);
p = struct();
p.s1 = struct();
p.s1.sp = subplot(1,2,1);
p.s1.p1 = struct();
p.s1.p2 = struct();
hold on
for i = 1:numel(E)
    ptmp = plot( nan*N,nan*E{i}(1,:),...
                        MarkerSize=sqrt(marker_size(i)),...
                        Marker=markers{i},...
                        MarkerEdgeColor=colors(i,:),...
                        Color=[colors(i,:),1],...
                        LineStyle=linspec{i});
    if (size(E{i},1)>1)
        N_ic = size(E{i},1);
        p.s1.p1(i).p = struct();
        p.s1.p2(i).p = struct();
        for j = 1:N_ic
            alpha = alpha_calc( j/N_ic );    
            p.s1.p1(i).p(j) = plot( N, E{i}(j,:),...
                                    Color=[colors(i,:),alpha],...
                                    LineStyle=linspec{i},...
                                    HandleVisibility='off');
            p.s1.p2(i).p(j) = scatter( N, E{i}(j,:),...
                                    marker_size(i),...
                                    markers{i},...
                                    MarkerEdgeColor=colors(i,:),...
                                    MarkerEdgeAlpha=alpha,...
                                    LineWidth=p.s1.p1(i).p(j).LineWidth,...
                                    HandleVisibility='off');
        end
    else
        p.s1.p1(i).p  = plot( N, E{i}(1,:),...
                              Color=[colors(i,:),1],...
                              LineStyle=linspec{i},...
                              HandleVisibility='off');
        p.s1.p2(i).p = scatter( N, E{i}(1,:),...
                                marker_size(i),...
                                markers{i},...
                                MarkerEdgeColor=colors(i,:),...
                                MarkerEdgeAlpha=1,...
                                LineWidth=p.s1.p1(i).p.LineWidth,...
                                HandleVisibility='off');
    end
    % plot(N,E{i},markers{i},MarkerSize=marker_size)
end
xlim([lim1,lim2])
xlabel(x_label1,Interpreter="latex");
ylabel(y_label1,Interpreter="latex");
p.s1.l = legend(legend_labels,Interpreter="latex");
p.s1.l.ItemTokenSize = [22,0];
set(gca,'Yscale','log')
set(gca,'Xscale','log')


p.s2 = struct();
p.s2.sp = subplot(1,2,2);
p.s2.p1 = struct();
p.s2.p2 = struct();
hold on
for i = 1:numel(OOA)
    ptmp = plot( nan*N,nan*OOA{i}(1,:),...
                        MarkerSize=sqrt(marker_size(i)),...
                        Marker=markers{i},...
                        MarkerEdgeColor=colors(i,:),...
                        Color=[colors(i,:),1],...
                        LineStyle=linspec{i});
    if (size(OOA{i},1)>1)
        N_ic = size(OOA{i},1);
        p.s2.p1(i).p = struct();
        p.s2.p2(i).p = struct();
        for j = 1:N_ic
            alpha = alpha_calc( j/N_ic );    
            p.s2.p1(i).p(j) = plot( N, OOA{i}(j,:),...
                                    Color=[colors(i,:),alpha],...
                                    LineStyle=linspec{i},...
                                    HandleVisibility='off');
            p.s2.p2(i).p(j) = scatter( N, OOA{i}(j,:),...
                                    marker_size(i),...
                                    markers{i},...
                                    MarkerEdgeColor=colors(i,:),...
                                    MarkerEdgeAlpha=alpha,...
                                    LineWidth=p.s2.p1(i).p(j).LineWidth,...
                                    HandleVisibility='off');
        end
    else
        p.s2.p1(i).p  = plot( N, OOA{i}(1,:),...
                              Color=[colors(i,:),1],...
                              LineStyle=linspec{i},...
                              HandleVisibility='off');
        p.s2.p2(i).p = scatter( N, OOA{i}(1,:),...
                                marker_size(i),...
                                markers{i},...
                                MarkerEdgeColor=colors(i,:),...
                                MarkerEdgeAlpha=1,...
                                LineWidth=p.s2.p1(i).p.LineWidth, ...
                                HandleVisibility='off');
    end
    % plot(N,OOA{i},markers{i},MarkerSize=marker_size)
end
xlim([lim1,lim2])
xlabel(x_label2,Interpreter="latex");
ylabel(y_label2,Interpreter="latex");
p.s2.l = legend(legend_labels,Interpreter="latex");
p.s2.l.ItemTokenSize = [22,0];
set(gca,'Xscale','log')
cellfun(@eval,post_plot_commands)
end

function alpha = alpha_calc(x)
min_alpha = 0.1;
max_alpha = 0.25;
alpha = min_alpha + (max_alpha-min_alpha)*sqrt(x);
% alpha = min_alpha + (1-min_alpha)*(2*x-1).^2;
end

function hfig = stdplot(i)
fontsize  = 6;%14;
linewidth = 1;%2;
% fontsize  = 20;
% linewidth = 2;
hfig=figure(i);
clf(hfig);
dim = [7.5 5.5 6.25 2.5];
set(hfig,'Units','Inches','Position',dim);
set(hfig,'DefaultAxesFontName','Helvetica');
set(hfig,'DefaultTextFontName','Helvetica'); 
set(hfig,'DefaultAxesFontSize',fontsize);
set(hfig,'DefaultTextFontSize',fontsize);
set(hfig,'PaperUnits',get(gcf,'Units'));
pos = get(hfig,'Position');
set(hfig,'PaperPosition',[0 0 pos(3) pos(4)]);
set(gca,'Units','Inches');
set(hfig,'DefaultLineLineWidth',linewidth)
set(hfig,'DefaultLineLineWidth',linewidth)

end