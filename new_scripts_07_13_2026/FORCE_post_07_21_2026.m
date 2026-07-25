%% Parsing KT-airfoil data (07/21/2026)
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
inputs.tau     = deg2rad(10.0);
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
foldernames1 = {};

% foldernames1 = [foldernames1,'JOUKOWSKI_C_GRID_curved_v2_04_08_regress_regress_2026-07-13_14.26.26'];
% foldernames1 = [foldernames1,'KT_AR_1_1_P2_ALPHA_05_2026-07-14_00.34.46'];
% foldernames1 = [foldernames1,'KT_AR_1_10_P2_ALPHA_05_2026-07-14_00.35.35'];
% foldernames1 = [foldernames1,'KT_AR_1_1_P4_ALPHA_05_2026-07-20_10.39.40'];
% foldernames1 = [foldernames1,'KT_AR_1_10_P4_ALPHA_05_2026-07-20_10.40.03'];
% foldernames1 = [foldernames1,'KT_AR_1_1_P4_ALPHA_00_2026-07-20_19.05.51'];
% foldernames1 = [foldernames1,'KT_AR_1_10_P4_ALPHA_00_2026-07-20_19.06.05'];
% foldernames1 = [foldernames1,'KT_AR_1_1_P4_ALPHA_05A_2026-07-20_20.45.05'];
% foldernames1 = [foldernames1,'KT_AR_1_1_P4_ALPHA_05A_N20_2026-07-20_22.23.14'];
% foldernames1 = [foldernames1,'KT_AR_1_1_TE_10_P4_ALPHA_05_2026-07-20_22.15.53'];
% foldernames1 = [foldernames1,'KT_AR_1_1_TE_10_P4_ALPHA_05_NEW_2026-07-21_13.39.19'];
% foldernames1 = [foldernames1,'KT_AR_1_1_P4_ALPHA_05A_D1_2026-07-21_16.31.38'];
% foldernames1 = [foldernames1,'KT_AR_1_1_P4_ALPHA_05A_D02_2026-07-22_14.15.01'];
% foldernames1 = [foldernames1,'KT_AR_1_1_P4_ALPHA_05A_D005_2026-07-22_14.15.18'];
% foldernames1 = [foldernames1,'KT_AR_1_1_P4_ALPHA_00_NEW_2026-07-23_10.23.07'];
% foldernames1 = [foldernames1,'KT_AR_1_1_P4_ALPHA_05_NEW_2026-07-23_20.13.17'];
% foldernames1 = [foldernames1,'KT_AR_1_1_TE_10_P4_ALPHA_05_NEW_2026-07-24_11.46.58'];
foldernames1 = [foldernames1,'KT_AR_1_1_TE_10_P4_ALPHA_05_NEW_PW_2026-07-24_18.46.06'];


foldernames1 = cellfun(@(str_b)strcat(DATA_DIR,str_b),foldernames1,UniformOutput=false);


folder = foldernames1{1};
S = get_soln_data_from_directory(folder);
G = get_grid_data_from_directory(folder);
DATA1 = get_airfoil_force_data_from_directory_alt(folder,inputs.alpha,inputs.nskip,airfoil,inputs.rho_ref,inputs.p_ref,inputs.a_ref,false);

% folder = foldernames1{2};
% S = get_soln_data_from_directory(folder);
% G = get_grid_data_from_directory(folder);
DATA2 = get_airfoil_force_data_from_directory_alt(folder,inputs.alpha,inputs.nskip,airfoil,inputs.rho_ref,inputs.p_ref,inputs.a_ref,true);

N1 = [DATA1.F(:).N];

% account for grid skip (find a better way to do this)
N1 = N1/inputs.nskip;

N_grids = numel(N1);

var = 'CL';
err_var_primal_1 = abs([DATA1.H(:).(['primal_',var])]-airfoil.(var));
err_var_primal_2 = abs([DATA2.H(:).(['primal_',var])]-airfoil.(var));
err_var_ete_1    = abs([DATA1.H(:).(['ete_',var])]-airfoil.(var));
err_var_ete_2    = abs([DATA2.H(:).(['ete_',var])]-airfoil.(var));

ooa_var_primal_1 = nan*err_var_primal_1;
ooa_var_primal_2 = nan*err_var_primal_2;
ooa_var_ete_1 = nan*err_var_ete_1;
ooa_var_ete_2 = nan*err_var_ete_2;

n_iter1 = size(ooa_var_ete_1,1);

r_fac = 2;
for i = 2:N_grids
    ooa_var_primal_1(:,i) = log(err_var_primal_1(:,i-1) ./ err_var_primal_1(:,i))./log(r_fac);
    ooa_var_primal_2(:,i) = log(err_var_primal_2(:,i-1) ./ err_var_primal_2(:,i))./log(r_fac);
    ooa_var_ete_1(:,i) = log(err_var_ete_1(:,i-1) ./ err_var_ete_1(:,i))./log(r_fac);
    ooa_var_ete_2(:,i) = log(err_var_ete_2(:,i-1) ./ err_var_ete_2(:,i))./log(r_fac);
end



%% Grid convergence of CL/CD error
print_ERR=false;
print_OOA=false;
% target_folder = 'C:\Users\wajordan\Desktop\CCAS_Annual_Review_Plots\OOA\CL';
target_folder = 'C:\Users\wajordan\Desktop\';
err_file = 'ERR_CL.png';
ooa_file = 'OOA_CL.png';

lim1 = 10^( floor( log10( N1(1)   ) ) );
lim2 = 10^( ceil(  log10( N1(end) ) ) );

hfig1 = stdplot(1);
subplot(1,2,1)
hold on

plot(N1,err_var_primal_1,'r-o',MarkerSize=3)
plot(N1,err_var_primal_2,'r--.',MarkerSize=3)
plot(N1,err_var_ete_1(1,:),'b-s',MarkerSize=3)
plot(N1,err_var_ete_1(end,:),'b:^',MarkerSize=3,HandleVisibility='off')
plot(N1,err_var_ete_2(1,:),'g--s',MarkerSize=3)
plot(N1,err_var_ete_2(end,:),'g-.^',MarkerSize=3,HandleVisibility='off')

xlim([lim1,lim2])

xlabel('$N_{cells}^{1/2}$','Interpreter','latex');
if strcmp(var,'CL')
    ylabel('$C_L$ Error','Interpreter','latex')
elseif strcmp(var,'CD')
    ylabel('$C_D$ Error','Interpreter','latex')
end
legend({'primal', 'primal (reconstructed)','corrected','corrected (reconstructed)'},Interpreter="latex")
set(gca,'Yscale','log')
set(gca,'Xscale','log')


subplot(1,2,2)
hold on
plot(N1,ooa_var_primal_1,'r-o',MarkerSize=3)
plot(N1,ooa_var_primal_2,'r--.',MarkerSize=3)
plot(N1,ooa_var_ete_1(1,:),'b-s',MarkerSize=3)
plot(N1,ooa_var_ete_1(end,:),'b:^',MarkerSize=3,HandleVisibility='off')
plot(N1,ooa_var_ete_2(1,:),'g--s',MarkerSize=3)
plot(N1,ooa_var_ete_2(end,:),'g-.^',MarkerSize=3,HandleVisibility='off')
xlim([lim1,lim2])
xlabel('$N_{cells}^{1/2}$','Interpreter','latex');
if strcmp(var,'CL')
    ylabel('$C_L$ OOA','Interpreter','latex')
elseif strcmp(var,'CD')
    ylabel('$C_D$ OOA','Interpreter','latex')
end
legend({'primal', 'primal (reconstructed)','corrected','corrected (reconstructed)'},Interpreter="latex")
set(gca,'Xscale','log')

post_plot_commands = {"set(hfig1.Children(4),'Ylim',[1e-6,1e0]);",...
                      "yticks(hfig1.Children(4),10.^(-6:1:0));",  ...
                      "set(hfig1.Children(4),'Xlim',[10,1000])",    ...
                      "set(hfig1.Children(3),'Location','southwest');",...
                      "set(hfig1.Children(2),'Ylim',[0,5]);",...
                      "set(hfig1.Children(2),'Xlim',[10,1000]);",...
                      "set(hfig1.Children(1),'Location','southwest');"};
cellfun(@eval,post_plot_commands)



if (print_ERR)
    exportgraphics(hfig1.Children(4),fullfile(target_folder,err_file),'Resolution',600)
end
if (print_OOA)
    exportgraphics(hfig1.Children(2),fullfile(target_folder,ooa_file),'Resolution',600)
end



ind = 6;

%% CL vs iterative corrections
var_ic1  = DATA1.H(ind).(['ic_',var]);

var_ic2  = DATA2.H(ind).(['ic_',var]);

n_iter = numel(var_ic1);

var_ex   = zeros(n_iter,1) + airfoil.(var);
var_ex1  = zeros(n_iter,1) + DATA1.H(ind).(['exact_',var]);
var_ex2  = zeros(n_iter,1) + DATA2.H(ind).(['exact_',var]);
var_pri1 = zeros(n_iter,1) + DATA1.H(ind).(['primal_',var]);
var_pri2 = zeros(n_iter,1) + DATA2.H(ind).(['primal_',var]);

figure(2);
hold on;
plot(var_ex,'k')
plot(var_ex1,'k--')
plot(var_ex2,'k:')
plot(var_pri1,'r')
plot(var_pri2,'r--')
plot(var_ic1, 'b--^')
plot(var_ic2, 'g--^')
legend({'analytic','exact (discrete)','exact (reconstructed)','primal','primal (reconstructed)','corrected','corrected (reconstructed)'},'Location','southeastoutside')
xlabel('iteration (iterative correction)')
if strcmp(var,'CL')
    ylabel('lift coefficient');
elseif strcmp(var,'CD')
    ylabel('drag coefficient');
end






%% figure 1: Cp Error (with respect to discrete exact Cp)
figure(3)
hold on
sz = numel(DATA1.F(ind).XC)/2;
top = 1:sz;
bot = sz+1:2*sz;

side = bot;
plot(  DATA1.F(ind).XC(side), ( DATA1.F(ind).primal_CP(side)     - DATA1.F(ind).exact_ana_CP(side) ),'r')
plot(  DATA2.F(ind).XC(side), ( DATA2.F(ind).primal_CP(side)     - DATA1.F(ind).exact_ana_CP(side) ),'r--')
plot(  DATA1.F(ind).XC(side), ( DATA1.F(ind).ic_CP(side,:)       - DATA1.F(ind).exact_ana_CP(side) ),'b')
plot(  DATA2.F(ind).XC(side), ( DATA2.F(ind).ic_CP(side,:)       - DATA2.F(ind).exact_ana_CP(side) ),'g')
% set(gca,'YScale','log')

legend({'primal','ETE1','ETE2'})
xlabel('x')
ylabel('Cp error')

function DE = calc_OOA(DE,r_fac)
N_grids = length(DE);
DE(1).OOA = nan*DE(1).E;
for i = 2:N_grids
    den = max(DE(i).E,eps(1));
    Etmp = DE(i-1).E ./ den;
    DE(i).OOA = log(Etmp)./log(r_fac);
end

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