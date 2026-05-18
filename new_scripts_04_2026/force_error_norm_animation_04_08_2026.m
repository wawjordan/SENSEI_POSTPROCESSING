%% Parsing KT-airfoil data (04/08/2026)
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

inputs = struct();
inputs.epsilon = 0.1;
inputs.kappa   = 0.0;
inputs.tau     = 0.0;
inputs.vinf    = 75.0;
inputs.rhoinf  = 1.0;
inputs.pinf    = 100000.0;
inputs.gamma   = 1.4;

inputs.alpha   = 5; % (degrees)
inputs.nskip   = 4;
inputs.rho_ref = 1.0;
inputs.p_ref   = 100000.0;
inputs.a_ref   = sqrt(inputs.gamma*inputs.p_ref/inputs.rho_ref);

% nondimensionalize inputs
inputs.vinf   = inputs.vinf/inputs.a_ref;
inputs.rhoinf = inputs.rhoinf/inputs.rho_ref;
inputs.pinf   = inputs.pinf/(inputs.rho_ref*inputs.a_ref^2);

airfoil        = kt_airfoil( inputs.epsilon, inputs.kappa, inputs.tau );
airfoil.vinf   = inputs.vinf;
airfoil.rhoinf = inputs.rhoinf;
airfoil.pinf   = inputs.pinf;
airfoil        = airfoil.set_alpha(inputs.alpha);

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

folder = foldernames1{5};
S = get_soln_data_from_directory(folder);
G = get_grid_data_from_directory(folder);
DATA1 = get_airfoil_force_data_from_directory_alt(folder,inputs.alpha,inputs.nskip,airfoil,inputs.rho_ref,inputs.p_ref,inputs.a_ref,false);
% DATA2 = get_airfoil_force_data_from_directory_alt(folder,inputs.alpha,inputs.nskip,airfoil,inputs.rho_ref,inputs.p_ref,inputs.a_ref,true);

folder = foldernames1{6};
S = get_soln_data_from_directory(folder);
G = get_grid_data_from_directory(folder);
DATA2 = get_airfoil_force_data_from_directory_alt(folder,inputs.alpha,inputs.nskip,airfoil,inputs.rho_ref,inputs.p_ref,inputs.a_ref,false);

N1 = [DATA1.F(:).N];
N2 = [DATA2.F(:).N];

% account for grid skip (find a better way to do this)
N1 = N1/inputs.nskip;

N_grids = numel(N1);

var = 'CD';
err_var_primal_1_1 = abs([DATA1.H(:).(['primal_',var])]-airfoil.(var));



%% Grid convergence of CL/CD error
print_ERR=true;
print_OOA=true;
target_folder = 'C:\Users\wajordan\Desktop\Meeting_Plots\ASME_VVUQ_2026_plots\ERR_CL';
v = VideoWriter(fullfile(target_folder,[var,'_error_anim_10.mp4']),"MPEG-4");
v.Quality = 100;
v.FrameRate = 2;
open(v);
new = true;
for j = 0:200
    % for i = 1:size(DATA1.H(i).(['ic_',var]))
    err_var_ete_1_1    = abs([DATA1.H(:).(['ic_',var])]-airfoil.(var));
    err_var_ete_2_1    = abs([DATA2.H(:).(['ic_',var])]-airfoil.(var));
    err_file = sprintf('ERR_CL_%0.2d.png',j);
    ooa_file = sprintf('OOA_CL_%0.2d.png',j);
    err_var_primal_1 = err_var_primal_1_1(1,:);
    err_var_ete_1    = err_var_ete_1_1(j+1,:);
    err_var_ete_2    = err_var_ete_2_1(j+1,:);
    
    ooa_var_primal_1 = nan*err_var_primal_1;
    ooa_var_ete_1 = nan*err_var_ete_1;
    ooa_var_ete_2 = nan*err_var_ete_2;
    
    r_fac = 2;
    for i = 2:N_grids
        ooa_var_primal_1(i) = log(err_var_primal_1(i-1) ./ err_var_primal_1(i))./log(r_fac);
        ooa_var_ete_1(i) = log(err_var_ete_1(i-1) ./ err_var_ete_1(i))./log(r_fac);
        ooa_var_ete_2(i) = log(err_var_ete_2(i-1) ./ err_var_ete_2(i))./log(r_fac);
    end
    
    lim1 = 10^( floor( log10( N1(1)   ) ) );
    lim2 = 10^( ceil(  log10( N1(end) ) ) );
    
    hfig1 = stdplot(1);
    subplot(1,2,1)
    hold on
    
    plot(N1,err_var_primal_1,'r-o',MarkerSize=3)
    plot(N1,err_var_ete_1,'b:s',MarkerSize=3)
    plot(N1,err_var_ete_2,'g--d',MarkerSize=3)
    
    xlim([lim1,lim2])
    
    xlabel('$N_{cells}^{1/2}$','Interpreter','latex');
    if strcmp(var,'CL')
        ylabel('$C_L$ Error','Interpreter','latex')
    elseif strcmp(var,'CD')
        ylabel('$C_D$ Error','Interpreter','latex')
    end
    legend({'Primal','ETE (K-exact)','ETE (CWENO)'},Interpreter="latex")
    set(gca,'Yscale','log')
    set(gca,'Xscale','log')
    
    
    subplot(1,2,2)
    hold on
    plot(N1,ooa_var_primal_1,'r-o',MarkerSize=3)
    plot(N1,ooa_var_ete_1,'b:s',MarkerSize=3)
    plot(N1,ooa_var_ete_2,'g--d',MarkerSize=3)
    xlim([lim1,lim2])
    xlabel('$N_{cells}^{1/2}$','Interpreter','latex');
    if strcmp(var,'CL')
        ylabel('$C_L$ OOA','Interpreter','latex')
    elseif strcmp(var,'CD')
        ylabel('$C_D$ OOA','Interpreter','latex')
    end
    legend({'Primal','ETE (K-exact)','ETE (CWENO)'},Interpreter="latex")
    set(gca,'Xscale','log')
    
    post_plot_commands = {"set(hfig1.Children(4),'Ylim',[1e-6,1e0]);",...
                          "yticks(hfig1.Children(4),10.^(-6:1:0));",  ...
                          "set(hfig1.Children(4),'Xlim',[10,1000])",    ...
                          "set(hfig1.Children(3),'Location','southwest');",...
                          "set(hfig1.Children(2),'Ylim',[0,5]);",...
                          "set(hfig1.Children(2),'Xlim',[10,1000]);",...
                          "set(hfig1.Children(1),'Location','southwest');"};
    cellfun(@eval,post_plot_commands)
    label_iteration(hfig1.Children(4),j);
    label_iteration(hfig1.Children(2),j);
    cellfun(@eval,post_plot_commands)
    F = print('-RGBImage','-r300');
    writeVideo(v,F)

    % extra frame in case matlab fucks it up
    if (j==0)
        writeVideo(v,F)
    end
    % if (print_ERR)
    %     % exportgraphics(hfig1.Children(4),fullfile(target_folder,err_file),'Resolution',600,Append=~new)
    %     exportgraphics(hfig1.Children(4),fullfile(target_folder,err_file),'Resolution',600)
    %     new = false;
    % end
    % if (print_OOA)
    %     % exportgraphics(hfig1.Children(2),fullfile(target_folder,ooa_file),'Resolution',600,Append=~new)
    %     exportgraphics(hfig1.Children(2),fullfile(target_folder,ooa_file),'Resolution',600)
    %     new = false;
    % end
end
% extra frame in case matlab fucks it up
writeVideo(v,F)
close(v);

function DE = calc_OOA(DE,r_fac)
N_grids = length(DE);
DE(1).OOA = nan*DE(1).E;
for i = 2:N_grids
    den = max(DE(i).E,eps(1));
    Etmp = DE(i-1).E ./ den;
    DE(i).OOA = log(Etmp)./log(r_fac);
end

end

function label_iteration(ax,iter)

xpos = 0.4;
ypos = 0.97;
txt = sprintf('Iterative Correction: %d',iter);

text(ax,xpos,ypos,txt,'Interpreter','latex','Units','Normalized',...
                      'VerticalAlignment','top')
end

function stdprint(hax,filename)
exportgraphics(hax,filename,'Resolution',600)
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