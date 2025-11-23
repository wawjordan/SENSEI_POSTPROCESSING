%% Force & Force History reading for grid convergence study prototype
% 02/12/2025
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

folder1 = 'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2025-02-18_08.52.56_OLD_STENCILS_UR=0.5_new_no_ho_exact_linear';
folder2 = 'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2025-02-17_15.51.22_OLD_STENCILS_UR=0.1_new';
% folder1 = 'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2025-02-14_08.43.08_OLD_STENCILS_limiters_on';
% folder1 = 'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2025-02-12_23.42.00_OLD_STENCILS_no_constraints';
% folder2 = 'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2025-02-13_10.12.07_OLD_STENCILS_rec=true';

% folder2 = 'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2025-02-17_17.29.45_OLD_STENCILS_UR=0.1_new_no_ho_exact';
% folder2 = 'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2_2025-02-12_09.01.45_NEW_STENCILS';
% folder2 = 'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2_2025-02-14_08.42.49_NEW_STENCILS_limiters_on';
% folder2 = 'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2_2025-02-12_23.41.00_NEW_STENCILS_no_constraints';
% folder2 = 'C:\Users\wajordan\Desktop\JOUKOWSKI_C_GRID_curved_2_2025-02-13_10.12.20_NEW_STENCILS_rec=true';

% S = get_soln_data_from_directory(folder2);
% G = get_grid_data_from_directory(folder2);

DATA1 = get_airfoil_force_data_from_directory_alt(folder1,alpha,nskip,airfoil);
DATA2 = get_airfoil_force_data_from_directory_alt(folder2,alpha,nskip,airfoil);
%%
% convert to Cp if needed
for i = 1:numel(DATA1.F)
    DATA1.F(i).primal_P = (DATA1.F(i).primal_P - pinf*(rho_ref*a_ref^2) ) ...
                        / (0.5*(rhoinf*rho_ref)*(vinf*a_ref)^2);
    DATA1.F(i).ete_P = (DATA1.F(i).ete_P - pinf*(rho_ref*a_ref^2) ) ...
                        / (0.5*(rhoinf*rho_ref)*(vinf*a_ref)^2);
    DATA1.F(i).ic_P = (DATA1.F(i).ic_P - pinf*(rho_ref*a_ref^2) ) ...
                        / (0.5*(rhoinf*rho_ref)*(vinf*a_ref)^2);
    DATA1.F(i).exact_P = (DATA1.F(i).exact_P - pinf*(rho_ref*a_ref^2) ) ...
                        / (0.5*(rhoinf*rho_ref)*(vinf*a_ref)^2);
    DATA1.F(i).exact_P2 = airfoil.get_averaged_cp_on_airfoil(DATA1.G(i).x,DATA1.G(i).y,true,true);
    DATA1.F(i).exact_P3 = airfoil.get_averaged_cp_on_airfoil(DATA1.G(i).x,DATA1.G(i).y,true,false);
end
for i = 1:numel(DATA2.F)
    DATA2.F(i).primal_P = (DATA2.F(i).primal_P - pinf*(rho_ref*a_ref^2) ) ...
                        / (0.5*(rhoinf*rho_ref)*(vinf*a_ref)^2);
    DATA2.F(i).ete_P = (DATA2.F(i).ete_P - pinf*(rho_ref*a_ref^2) ) ...
                        / (0.5*(rhoinf*rho_ref)*(vinf*a_ref)^2);
    DATA2.F(i).ic_P = (DATA2.F(i).ic_P - pinf*(rho_ref*a_ref^2) ) ...
                        / (0.5*(rhoinf*rho_ref)*(vinf*a_ref)^2);
    DATA2.F(i).exact_P = (DATA2.F(i).exact_P - pinf*(rho_ref*a_ref^2) ) ...
                        / (0.5*(rhoinf*rho_ref)*(vinf*a_ref)^2);
    DATA2.F(i).exact_P2 = airfoil.get_averaged_cp_on_airfoil(DATA2.G(i).x,DATA2.G(i).y,true,true);
    DATA2.F(i).exact_P3 = airfoil.get_averaged_cp_on_airfoil(DATA2.G(i).x,DATA2.G(i).y,true,false);
end

err_CL_primal_1 = abs([DATA1.H(:).primal_CL]-airfoil.CL);

err_CL_ete_1 = abs([DATA1.H(:).ete_CL]-airfoil.CL);
err_CL_ete_2 = abs([DATA2.H(:).ete_CL]-airfoil.CL);
N1 = [DATA1.F(:).N];
N2 = [DATA2.F(:).N];



%% Grid convergence of CL error
figure(1)
hold on
plot(N1,err_CL_primal_1(1,:),'r-o')
plot(N1,err_CL_ete_1(1,:),'b-^')
plot(N1,err_CL_ete_1(2,:),'b--^')
plot(N2,err_CL_ete_2(1,:),'g-v')
plot(N2,err_CL_ete_2(2,:),'g--v')
set(gca,'Yscale','log')
set(gca,'Xscale','log')
legend({'primal','ETE (old stencil)','IC 150 (old stencil)','ETE (new stencil)','IC 150 (new stencil)'})



ind = 5;
%% figure 2: Cp plots
figure(2)
hold on
[h1,h2] = airfoil.plot_cp;
h1.Color='k';
h2.Color='k';
h1.HandleVisibility='off';
plot(DATA1.F(ind).XC,DATA1.F(ind).exact_P,'k--')
plot(DATA1.F(ind).XC,DATA1.F(ind).exact_P2,'m-')
plot(DATA1.F(ind).XC,DATA1.F(ind).exact_P3,'m--')
plot(DATA1.F(ind).XC,DATA1.F(ind).primal_P,'r')
plot(DATA1.F(ind).XC,DATA1.F(ind).ete_P(:,1),'b')
plot(DATA1.F(ind).XC,DATA1.F(ind).ete_P(:,2),'b--')
plot(DATA2.F(ind).XC,DATA2.F(ind).ete_P(:,1),'g')
plot(DATA2.F(ind).XC,DATA2.F(ind).ete_P(:,2),'g--')
legend({'analytic','exact (discrete)','exact (linear boundary)','exact (true boundary)','primal','ETE (old stencil)','IC 150 (old stencil)','ETE (new stencil)','IC 150 (new stencil)'})
xlabel('x/c')
ylabel('Cp')
%% figure 3: Cp Error (with respect to discrete exact Cp)
figure(3)
hold on
plot(DATA1.F(ind).primal_P   - DATA1.F(ind).exact_P2,'r')
plot(DATA1.F(ind).ete_P(:,1) - DATA1.F(ind).exact_P2,'g')
plot(DATA1.F(ind).ete_P(:,2) - DATA1.F(ind).exact_P2,'g--')
plot(DATA2.F(ind).ete_P(:,1) - DATA2.F(ind).exact_P2,'b')
plot(DATA2.F(ind).ete_P(:,2) - DATA2.F(ind).exact_P2,'b--')
legend({'primal','ETE (old stencil)','IC=150 (old stencil)','ETE (new stencil)','IC=150 (new stencil)'})
xlabel('i index')
ylabel('Cp error')

plot(DATA1.S(ind).src_err_lin,'m')
plot(DATA1.S(ind).src_err,'m--')
%% figure 4: Cp absolute Error (with respect to discrete exact Cp)
figure(4)
hold on
plot(abs( DATA1.F(ind).primal_P   - DATA1.F(ind).exact_P2 ),'r')
plot(abs( DATA1.F(ind).ete_P(:,1) - DATA1.F(ind).exact_P2 ),'g')
plot(abs( DATA1.F(ind).ete_P(:,2) - DATA1.F(ind).exact_P2 ),'g--')
plot(abs( DATA2.F(ind).ete_P(:,1) - DATA2.F(ind).exact_P2 ),'b')
plot(abs( DATA2.F(ind).ete_P(:,2) - DATA2.F(ind).exact_P2 ),'b--')

plot(abs(DATA1.S(ind).src_err_lin),'m')
plot(abs(DATA1.S(ind).src_err),'m--')
set(gca,'YScale','log')
legend({'primal','ETE (old stencil)','IC=150 (old stencil)','ETE (new stencil)','IC=150 (new stencil)'})
xlabel('i index')
ylabel('|Cp| error')
%% figure 5: CL vs iterative corrections
CL_ic1  = DATA1.H(ind).ic_CL;
CL_ic2  = DATA2.H(ind).ic_CL;
n_iter = max(numel(CL_ic1),numel(CL_ic2));

CL_ex   = zeros(n_iter,1) + airfoil.CL;
CL_ex1  = zeros(n_iter,1) + DATA1.H(ind).exact_CL;
CL_ex2  = zeros(n_iter,1) + DATA2.H(ind).exact_CL;
CL_pri1 = zeros(n_iter,1) + DATA1.H(ind).primal_CL;
CL_pri2 = zeros(n_iter,1) + DATA2.H(ind).primal_CL;
CL_ete1 = zeros(n_iter,1) + DATA1.H(ind).ete_CL(1);
CL_ete2 = zeros(n_iter,1) + DATA2.H(ind).ete_CL(1);

figure(5);
hold on;
plot(CL_ex,'k')
plot(CL_ex1,'k--')
% plot(CL_ex2,'k--')
plot(CL_pri1,'r')
plot(CL_ete1, 'b')
plot(CL_ic1, 'b--^')
% plot(CL_pri2,'r.-')
plot(CL_ete2, 'g')
plot(CL_ic2, 'g--v')
legend({'analytic','exact (discrete)','primal','ETE (old stencil)','IC (old stencil)','ETE (new stencil)','IC (new stencil)'},'Location','southeastoutside')
xlabel('iteration (iterative correction)')
ylabel('lift coefficient')
% %%
% CX_ic1  = DATA1.H(ind).ic_CX;
% CX_ic2  = DATA2.H(ind).ic_CX;
% n_iter = max(numel(CX_ic1),numel(CX_ic2));
% 
% CX_ex   = zeros(n_iter,1) + airfoil.CX;
% CX_ex1  = zeros(n_iter,1) + DATA1.H(ind).exact_CX;
% CX_ex2  = zeros(n_iter,1) + DATA2.H(ind).exact_CX;
% CX_pri1 = zeros(n_iter,1) + DATA1.H(ind).primal_CX;
% CX_pri2 = zeros(n_iter,1) + DATA2.H(ind).primal_CX;
% CX_ete1 = zeros(n_iter,1) + DATA1.H(ind).ete_CX(1);
% CX_ete2 = zeros(n_iter,1) + DATA2.H(ind).ete_CX(1);
% clf;
% hold on;
% plot(CX_ex,'k')
% plot(CX_ex1,'k--')
% plot(CX_ex2,'k--')
% plot(CX_pri1,'r')
% plot(CX_ete1, 'b--')
% plot(CX_ic1, 'b')
% plot(CX_pri2,'r.-')
% plot(CX_ete2, 'b:')
% plot(CX_ic2, 'b.-')
% 
% %%
% CY_ic1  = DATA1.H(ind).ic_CY;
% CY_ic2  = DATA2.H(ind).ic_CY;
% n_iter = max(numel(CY_ic1),numel(CY_ic2));
% 
% CY_ex   = zeros(n_iter,1) + airfoil.CY;
% CY_ex1  = zeros(n_iter,1) + DATA1.H(ind).exact_CY;
% CY_ex2  = zeros(n_iter,1) + DATA2.H(ind).exact_CY;
% CY_pri1 = zeros(n_iter,1) + DATA1.H(ind).primal_CY;
% CY_pri2 = zeros(n_iter,1) + DATA2.H(ind).primal_CY;
% CY_ete1 = zeros(n_iter,1) + DATA1.H(ind).ete_CY(1);
% CY_ete2 = zeros(n_iter,1) + DATA2.H(ind).ete_CY(1);
% clf;
% hold on;
% plot(CY_ex,'k')
% plot(CY_ex1,'k--')
% plot(CY_ex2,'k--')
% plot(CY_pri1,'r')
% plot(CY_ete1, 'b--')
% plot(CY_ic1, 'b')
% plot(CY_pri2,'r.-')
% plot(CY_ete2, 'b:')
% plot(CY_ic2, 'b.-')