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

inputs = struct();
inputs.epsilon = 0.1;
inputs.kappa   = 0.0;
inputs.tau     = 0.0;
inputs.vinf    = 75.0;
inputs.rhoinf  = 1.0;
inputs.pinf    = 100000.0;
inputs.gamma   = 1.4;

inputs.alpha   = 5; % (degrees)
inputs.nskip   = 2;
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

DATA_DIR='C:\Users\wajordan\Desktop\';
foldernames1 = {             'ALPHA_0_JOUKOWSKI_C_GRID_curved_2025-11-23_13.28.07_K_EXACT'};% 1
foldernames1 = [foldernames1,'ALPHA_0_JOUKOWSKI_C_GRID_curved_2025-11-23_12.47.49_VAR_REC_8_100'];% 2
foldernames1 = [foldernames1,'ALPHA_0_JOUKOWSKI_C_GRID_curved_2025-11-23_14.30.04_VAR_REC_8_1000'];% 3
foldernames1 = [foldernames1,'ALPHA_0_JOUKOWSKI_C_GRID_curved_2025-11-23_15.35.53_VAR_REC_ALL_100']; % 4
foldernames1 = [foldernames1,'ALPHA_0_JOUKOWSKI_C_GRID_curved_2025-11-23_15.03.36_VAR_REC_ALL_1000'];% 5
foldernames1 = [foldernames1,'ALPHA_0_JOUKOWSKI_C_GRID_curved_2025-11-23_15.57.49_VAR_REC_ALL_100_USE_HO_GRID_RECONSTRUCT=T'];% 6
foldernames1 = [foldernames1,'ALPHA_0_JOUKOWSKI_C_GRID_curved_2025-11-23_16.13.48_VAR_REC_ALL_100_no_bc'];% 7
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-11-26_15.23.56_K_EXACT_no_bc'];% 8
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-11-26_15.40.39_K_EXACT'];% 9
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-11-23_16.46.56_VAR_REC_ALL_100_no_bc'];% 10
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-11-26_16.20.39_VAR_REC_ALL_1000_no_bc'];% 11
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-11-26_15.54.04_VAR_REC_ALL_100_bc'];% 12

foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-04_12.39.32_KEXACT_yes_bc_IC_10_ur0-5'];% 13


foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-01_12.50.36_VAR_REC_ALL_1000_USE_HO_GRID_RECONSTRUCT=T_no_bc_IC_10_no_ur'];% 14
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-01_23.04.15_VAR_REC_ALL_1000_no_bc_IC_10_ur0-5'];% 15
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-02_10.32.10_VAR_REC_ALL_10_no_bc_IC_10_ur0-5'];% 16

foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-02_19.07.11_VAR_REC_1_1000_no_bc_IC_10_ur0-5'];% 17
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-02_16.26.38_VAR_REC_2_1000_no_bc_IC_10_ur0-5'];% 18
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-02_14.01.27_VAR_REC_4_1000_no_bc_IC_10_ur0-5'];% 19
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-02_12.18.32_VAR_REC_8_1000_no_bc_IC_10_ur0-5'];% 20

foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-03_20.39.41_VAR_REC_8_10_yes_bc_IC_10_ur0-5'];% 21
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-04_11.01.38_VAR_REC_8_100_yes_bc_IC_10_ur0-5'];% 22
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-03_11.55.12_VAR_REC_8_1000_yes_bc_IC_10_ur0-5'];% 23
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-05_11.19.30_VAR_REC_1_100_yes_bc_IC_10_ur0-5'];% 24
% abe23e21f0a6d7f1466407415ad2ebb6f5925edc


foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-02_17.52.48_VAR_REC_8_100_no_bc_IC_200_ur0-5'];% 25
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-08_13.39.39_VAR_REC_ALL_1000_no_bc_IC_200_ur0-5'];% 26


%% New compile (12d0c99f7af602e3eac64ad5bde538f3262efefb)
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-05_13.20.34_VAR_REC_1_100_yes_bc_IC_10_ur0-5_NEW'];% 27
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-05_18.06.19_VAR_REC_1_100_yes_bc_IC_10_ur0-5_NEW_10_iter'];% 28
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-08_13.32.09_VAR_REC_2_100_yes_bc_IC_10_ur0-5_NEW_10_iter'];% 29
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-09_13.43.19_ORDER_3_VAR_REC_1_100_yes_bc_IC_10_ur0-5_NEW_10_iter'];% 30
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-09_19.44.01_ORDER_3_VAR_REC_1_100_no_bc_IC_10_ur0-5_NEW_10_iter'];% 31
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_newer_2025-12-09_21.41.18_ORDER_3_VAR_REC_1_100_no_bc_IC_10_ur0-5_NEWER_10_iter'];% 32
foldernames1 = [foldernames1,'ALPHA_5_JOUKOWSKI_C_GRID_curved_2025-12-10_10.50.20_ORDER_3_VAR_REC_1_100_no_bc_IC_10_ur0-5_NEW_10_iter_GEO4'];% 33

foldernames1 = cellfun(@(str_b)strcat(DATA_DIR,str_b),foldernames1,UniformOutput=false);

folder = foldernames1{17};
S = get_soln_data_from_directory(folder);
G = get_grid_data_from_directory(folder);
DATA1 = get_airfoil_force_data_from_directory_alt(folder,inputs.alpha,inputs.nskip,airfoil,inputs.rho_ref,inputs.p_ref,inputs.a_ref,false);
% DATA2 = get_airfoil_force_data_from_directory_alt(folder,inputs.alpha,inputs.nskip,airfoil,inputs.rho_ref,inputs.p_ref,inputs.a_ref,true);

folder = foldernames1{20};
S = get_soln_data_from_directory(folder);
G = get_grid_data_from_directory(folder);
DATA2 = get_airfoil_force_data_from_directory_alt(folder,inputs.alpha,inputs.nskip,airfoil,inputs.rho_ref,inputs.p_ref,inputs.a_ref,false);

N1 = [DATA1.F(:).N];
N2 = [DATA2.F(:).N];

% account for grid skip (find a better way to do this)
N1 = N1/2;

N_grids = numel(N1);

var = 'CL';
err_var_primal_1 = abs([DATA1.H(:).(['primal_',var])]-airfoil.(var));
err_var_ete_1    = abs([DATA1.H(:).(['ete_',var])]-airfoil.(var));
err_var_ete_2    = abs([DATA2.H(:).(['ete_',var])]-airfoil.(var));

% don't worry about ic for the moment
err_var_primal_1 = err_var_primal_1(1,:);
err_var_ete_1    = err_var_ete_1(1,:);
err_var_ete_2    = err_var_ete_2(1,:);

ooa_var_primal_1 = nan*err_var_primal_1;
ooa_var_ete_1 = nan*err_var_ete_1;
ooa_var_ete_2 = nan*err_var_ete_2;

r_fac = 2;
for i = 2:N_grids
    ooa_var_primal_1(i) = log(err_var_primal_1(i-1) ./ err_var_primal_1(i))./log(r_fac);
    ooa_var_ete_1(i) = log(err_var_ete_1(i-1) ./ err_var_ete_1(i))./log(r_fac);
    ooa_var_ete_2(i) = log(err_var_ete_2(i-1) ./ err_var_ete_2(i))./log(r_fac);
end



%% Grid convergence of CL/CD error
print_ERR=false;
print_OOA=false;
target_folder = 'C:\Users\wajordan\Desktop\CCAS_Annual_Review_Plots\OOA\CL';
err_file = 'ERR.png';
ooa_file = 'OOA.png';

% figure(2)
% hold on
% plot(N1,err_var_primal_1(1,:),'r-s')
% plot(N1,err_var_ete_1(1,:),'b-^')
% for i = 1:size(err_var_ete_1,1)
%     plot(N1,err_var_ete_1(i,:),'r-')
% end
% plot(N1,err_var_ete_1(end,:),'b--^')
% plot(N1,err_var_ete_2(1,:),'g-s')
% for i = 1:size(err_var_ete_2,1)
%     plot(N1,err_var_ete_2(i,:),'b-')
% end
% plot(N1,err_var_ete_2(end,:),'g--s')
% set(gca,'Yscale','log')
% set(gca,'Xscale','log')
% legend({'primal','ETE (Old Rec.)','ETE (New Rec.)'})
% xlabel('')
% if strcmp(var,'CL')
%     ylabel('lift coefficient');
% elseif strcmp(var,'CD')
%     ylabel('drag coefficient');
% end
% markers = {'o','s','d','^'};


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
legend({'Primal','ETE (Old Rec.)','ETE (New Rec.)'},Interpreter="latex")
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
legend({'Primal','ETE (Old Rec.)','ETE (New Rec.)'},Interpreter="latex")
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


ind = 5;

%% figure 1: Cp Error (with respect to discrete exact Cp)
figure(1)
hold on
% plot(DATA1.F(ind).exact_lin_CP   - DATA1.F(ind).exact_ana_CN1P,'k')
% plot(DATA1.F(ind).exact_sim_CP   - DATA1.F(ind).exact_ana_CP,'k--')
% plot(DATA1.F(ind).primal_CP      - DATA1.F(ind).exact_ana_CP,'r')
% plot(DATA1.F(ind).ete_CP(:,1)    - DATA1.F(ind).exact_ana_CP,'g')
% plot(DATA1.F(ind).ete_CP(:,end)  - DATA1.F(ind).exact_ana_CP,'m')
% 
% plot(DATA2.F(ind).ete_CP(:,1)    - DATA2.F(ind).exact_ana_CP,'g--')
% plot(DATA2.F(ind).ete_CP(:,end)  - DATA2.F(ind).exact_ana_CP,'m--')
% 
% legend({'ex\_lin','ex\_sim','primal','ETE'})
% xlabel('i index')
% ylabel('Cp error')
sz = numel(DATA1.F(ind).XC)/2;
top = 1:sz;
bot = sz+1:2*sz;

side = top;
% plot(  DATA1.F(ind).XC(side), DATA1.F(ind).exact_lin_CP(side)   - DATA1.F(ind).exact_ana_CP(side),'k')
% plot(  DATA1.F(ind).XC(side), DATA1.F(ind).exact_sim_CP(side)   - DATA1.F(ind).exact_ana_CP(side),'k--')
plot(  DATA1.F(ind).XC(side), ( DATA1.F(ind).primal_CP(side)      - DATA1.F(ind).exact_ana_CP(side) ),'r')
plot(  DATA1.F(ind).XC(side), ( DATA1.F(ind).ete_CP(side,1)    - DATA1.F(ind).exact_ana_CP(side) ),'g')
% plot(  DATA1.F(ind).XC(side), ( DATA1.F(ind).ete_CP(side,end)  - DATA1.F(ind).exact_ana_CP(side) ),'m')

plot(  DATA2.F(ind).XC(side), ( DATA2.F(ind).ete_CP(side,1)    - DATA2.F(ind).exact_ana_CP(side) ),'g--')
% plot(  DATA2.F(ind).XC(side), ( DATA2.F(ind).ete_CP(side,end)  - DATA2.F(ind).exact_ana_CP(side) ),'m--')
% set(gca,'YScale','log')

legend({'primal','ETE (Old Rec.)','ETE (New Rec.)'})
xlabel('x')
ylabel('Cp error')


%% CL vs iterative corrections
% var_ic1  = DATA1.H(ind).(['ic_',var]);
% 
% var_ic2  = DATA2.H(ind).(['ic_',var]);
% 
% n_iter = numel(var_ic1);
% 
% var_ex   = zeros(n_iter,1) + airfoil.(var);
% var_ex1  = zeros(n_iter,1) + DATA1.H(ind).(['exact_',var]);
% var_ex2  = zeros(n_iter,1) + DATA2.H(ind).(['exact_',var]);
% var_pri1 = zeros(n_iter,1) + DATA1.H(ind).(['primal_',var]);
% 
% figure(3);
% hold on;
% plot(var_ex,'k')
% plot(var_ex1,'k--')
% plot(var_ex2,'k:')
% plot(var_pri1,'r')
% plot(var_ic1, 'b--^')
% plot(var_ic2, 'g--^')
% legend({'analytic','exact (discrete)','exact (discrete rec)','primal','IC'},'Location','southeastoutside')
% xlabel('iteration (iterative correction)')
% if strcmp(var,'CL')
%     ylabel('lift coefficient');
% elseif strcmp(var,'CD')
%     ylabel('drag coefficient');
% end


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