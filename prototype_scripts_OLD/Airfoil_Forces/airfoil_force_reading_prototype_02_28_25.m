%% Force & Force History reading for grid convergence study prototype
% 02/28/2025
clc; clear; close all; fclose('all');

epsilon = 0.1;
kappa   = 0.0;
tau     = 0.0;
vinf    = 75.0;
rhoinf  = 1.0;
pinf    = 100000.0;
gamma   = 1.4;

alpha   = 5; % (degrees)
nskip   = 2;
rho_ref = 1.0;
p_ref   = 100000.0;
a_ref   = sqrt(gamma*p_ref/rho_ref);

% nondimensionalize inputs
vinf   = vinf/a_ref;
rhoinf = rhoinf/rho_ref;
pinf   = pinf/(rho_ref*a_ref^2);

airfoil        = kt_airfoil(epsilon,kappa,tau);
airfoil.vinf   = vinf;
airfoil.rhoinf = rhoinf;
airfoil.pinf   = pinf;
airfoil        = airfoil.set_alpha(alpha);

% folder1 = 'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2025-02-28_09.01.38';
% folder1 = 'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2025-02-28_11.01.29_unconstrained';
% folder2 = 'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2025-03-05_16.57.02';
folder2 = 'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2_2025-11-13_12.06.39';
S = get_soln_data_from_directory(folder2);
G = get_grid_data_from_directory(folder2);

DATA1 = get_airfoil_force_data_from_directory_alt(folder2,alpha,nskip,airfoil,rho_ref,p_ref,a_ref,false);
DATA2 = get_airfoil_force_data_from_directory_alt(folder2,alpha,nskip,airfoil,rho_ref,p_ref,a_ref,true);

var = 'CL';
err_CL_primal_1 = abs([DATA1.H(:).(['primal_',var])]-airfoil.(var));
err_CL_primal_2 = abs([DATA2.H(:).(['primal_',var])]-airfoil.(var));

err_CL_ete_1 = abs([DATA1.H(:).(['ete_',var])]-airfoil.(var));
err_CL_ete_2 = abs([DATA2.H(:).(['ete_',var])]-airfoil.(var));
N1 = [DATA1.F(:).N];
N2 = [DATA2.F(:).N];

ind = 3;

% fprintf('Size = %i')

%% figure 2: Cp
% figure(2)
% hold on
% plot(DATA1.F(ind).primal_CP   - DATA1.F(ind).exact_ana_CP,'r')
% plot(DATA1.F(ind).ete_CP(:,1) - DATA1.F(ind).exact_ana_CP,'g')
% plot(DATA1.F(ind).ete_CP(:,2) - DATA1.F(ind).exact_ana_CP,'g--')
% plot(DATA2.F(ind).primal_CP   - DATA2.F(ind).exact_ana_CP,'r--')
% plot(DATA2.F(ind).ete_CP(:,1) - DATA2.F(ind).exact_ana_CP,'b')
% plot(DATA2.F(ind).ete_CP(:,2) - DATA2.F(ind).exact_ana_CP,'b--')
% legend({'primal','ETE','IC 50','primal (rec)','ETE (rec)','IC 50 (rec)'})
% xlabel('i index')
% ylabel('Cp error')

%% figure 3: Cp Error (with respect to discrete exact Cp)
figure(3)
hold on
plot(DATA1.F(ind).primal_CP   - DATA1.F(ind).exact_ana_CP,'r')
plot(DATA1.F(ind).ete_CP(:,1) - DATA1.F(ind).exact_ana_CP,'g')
plot(DATA1.F(ind).ete_CP(:,2) - DATA1.F(ind).exact_ana_CP,'g--')
plot(DATA2.F(ind).primal_CP   - DATA2.F(ind).exact_ana_CP,'r--')
plot(DATA2.F(ind).ete_CP(:,1) - DATA2.F(ind).exact_ana_CP,'b')
plot(DATA2.F(ind).ete_CP(:,2) - DATA2.F(ind).exact_ana_CP,'b--')
legend({'primal','ETE','IC 50','primal (rec)','ETE (rec)','IC 50 (rec)'})
xlabel('i index')
ylabel('Cp error')

% plot(DATA1.S(ind).src_err_lin,'m')
% plot(DATA1.S(ind).src_err,'m--')

%% Grid convergence of CL error
figure(1)
hold on
plot(N1,err_CL_primal_1(1,:),'r-s')
plot(N1,err_CL_ete_1(1,:),'b-^')
plot(N1,err_CL_ete_1(2,:),'b--^')

plot(N2,err_CL_primal_2(1,:),'r-o')
plot(N2,err_CL_ete_2(1,:),'g-v')
plot(N2,err_CL_ete_2(2,:),'g--v')

set(gca,'Yscale','log')
set(gca,'Xscale','log')
legend({'primal','ETE','IC 50','primal (rec)','ETE (rec)','IC 50 (rec)'})

%% CL vs iterative corrections
CL_ic1  = DATA1.H(ind).(['ic_',var]);
CL_ic2  = DATA2.H(ind).(['ic_',var]);
n_iter = max(numel(CL_ic1),numel(CL_ic2));

CL_ex   = zeros(n_iter,1) + airfoil.(var);
CL_ex1  = zeros(n_iter,1) + DATA1.H(ind).(['exact_',var]);
CL_ex2  = zeros(n_iter,1) + DATA2.H(ind).(['exact_',var]);
CL_pri1 = zeros(n_iter,1) + DATA1.H(ind).(['primal_',var]);
CL_pri2 = zeros(n_iter,1) + DATA2.H(ind).(['primal_',var]);

figure(2);
hold on;
plot(CL_ex,'k')
plot(CL_ex1,'k--')
plot(CL_ex2,'k:')
plot(CL_pri1,'r')
plot(CL_pri2,'r.-')
plot(CL_ic1, 'b--^')
plot(CL_ic2, 'g--v')

legend({'analytic','exact (discrete)','exact (discrete rec)','primal','primal rec','IC','IC (rec)'},'Location','southeastoutside')
xlabel('iteration (iterative correction)')
ylabel('lift coefficient')