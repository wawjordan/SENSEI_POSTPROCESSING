function DATA = get_airfoil_force_data_from_directory_alt(folder,alpha,nskip,airfoil,rho_ref,p_ref,a_ref,reconstruct)


if ( nargin > 7 )
    DATA.F = read_surface_forces_from_folder(folder,reconstruct);
    DATA.H = read_force_histories_from_folder(folder,alpha,reconstruct);
else
    DATA.F = read_surface_forces_from_folder(folder);
    DATA.H = read_force_histories_from_folder(folder,alpha);
end
DATA.G = read_surface_from_folder(folder,nskip);
DATA.S = get_source_terms_from_folder(folder,airfoil);

% calculate exact pressure integrals

p_off   = p_ref;
p_scale = (0.5*(airfoil.rhoinf*rho_ref)*(airfoil.vinf*a_ref)^2);

for i = 1:numel(DATA.F)
    DATA.F(i).primal_CP = ( DATA.F(i).primal_P - p_off ) / p_scale;
    DATA.F(i).ete_CP    = ( DATA.F(i).ete_P    - p_off ) / p_scale;
    DATA.F(i).ic_CP     = ( DATA.F(i).ic_P     - p_off ) / p_scale;
    DATA.F(i).exact_sim_CP  = ( DATA.F(i).exact_P  - p_off ) / p_scale;
    DATA.F(i).exact_lin_CP = airfoil.get_averaged_cp_on_airfoil(DATA.G(i).x,DATA.G(i).y,true,true);
    DATA.F(i).exact_ana_CP = airfoil.get_averaged_cp_on_airfoil(DATA.G(i).x,DATA.G(i).y,true,false);
end
end

function S = get_source_terms_from_folder(folder,airfoil)
D = get_soln_data_from_directory(folder);
S = struct();
N_files = length(D);
% get the corresponding source term data
tol = 100*eps(1);
for n = 1:N_files
    S(n).N  = D(n).N;
    x   = D(n).DATA.ZONE.X;
    y   = D(n).DATA.ZONE.Y;
    src = D(n).DATA.ZONE.SRC_ENERGY;
    mask = x(:,1,1)<=1+tol;
    idx = find(mask);
    S(n).src = src(idx(1:end-1),1,1);

    x1 = x(idx,1:2,1);
    y1 = y(idx,1:2,1);


    S(n).src_ex_lin = S(n).src * 0;
    S(n).src_ex     = S(n).src * 0;
    S(n).src_err_lin = S(n).src * 0;
    S(n).src_err     = S(n).src * 0;
    for i = 1:sum(mask)-1
        xtmp = [x1(i,1),x1(i+1,1),x1(i+1,2),x1(i,2)];
        ytmp = [y1(i,1),y1(i+1,1),y1(i+1,2),y1(i,2)];
        src_tmp = airfoil.calculate_source(xtmp,ytmp,true,false);
        S(n).src_ex_lin(i)  = src_tmp(5);
        src_tmp = airfoil.calculate_source(xtmp,ytmp,true,true);
        S(n).src_ex(i)  = src_tmp(5);
        S(n).src_err_lin(i)  = S(n).src(i) - S(n).src_ex_lin(i);
        S(n).src_err(i)  = S(n).src(i) - S(n).src_ex(i);
    end
end

end

function G = read_surface_from_folder(folder,nskip)

G1 = get_grid_data_from_directory(folder);
G = struct();
N_files = length(G1);
tol = 100*eps(1);
% get the corresponding source term data
for i = 1:N_files
    G(i).N  = G1(i).N;
    if G1(i).GRID.nblocks>1
        error('this function is only set up for single block grids')
    end
    x   = G1(i).GRID.gblock(1).x;
    y   = G1(i).GRID.gblock(1).y;
    mask = x(:,1,1)<=1+tol;
    G(i).x0 = x(mask,1,1);
    G(i).y0 = y(mask,1,1);
    G(i).x = G(i).x0(1:nskip:end);
    G(i).y = G(i).y0(1:nskip:end);
    G(i).xc = 0.5*( G(i).x(1:end-1) + G(i).x(2:end) );
    G(i).yc = 0.5*( G(i).y(1:end-1) + G(i).y(2:end) );
end
end

function F = read_surface_forces_from_folder(folder,reconstruct)
if ( nargin > 1 )&&( reconstruct )
    surface_pat = { '*-primal_rec-surface_forces.dat',    ...%
                    '*-ete_rec-surface_forces.dat',       ...%
                    '*-ic_rec-surface_forces.dat',        ...%
                    '*-exact_rec-surface_forces.dat' };
else
    surface_pat = { '*-primal-surface_forces.dat',    ...%
                   '*-ete-surface_forces.dat',       ...%
                   '*-ic-surface_forces.dat',        ...%
                   '*-exact-surface_forces.dat' };
end

F = struct();

% first check for the number of primal-force_history files
% (these will still be output, even if the iterative correction fails) 
TMP = dir(fullfile(folder,'**\',surface_pat{1}));
% for now, we will assume a single file per folder
N_files = length(TMP);
if (N_files==0)
    error('no matching primal files in folder')
end

folders = cell(N_files,1);
% get the corresponding folders and primal force data
for i = 1:N_files
    folders{i} = replace(TMP(i).folder,{folder,surface_pat{1}},{'',''});
    % fprintf('found folder: %s\n',folders{i});
    F(i).N  = cellfun(@str2double,regexp(folders{i},'\d*','match'));
    F(i).N  = sqrt(prod(F(i).N));
    file_name = fullfile(TMP(i).folder,TMP(i).name);
    [F(i).XC,F(i).primal_P] = get_airfoil_surface_pressure_data_from_file(file_name);
end

% Now loop through the folders and grab the corresponding files

for i = 1:N_files
    F(i).ete_P = nan*F(i).primal_P;
    TMP = dir(fullfile(folder,folders{i},surface_pat{2}));
    if (isempty(TMP))
        warning('No ete-surface_forces.dat in folder: %s\n',folders{i})
        continue
    end
    file_name = fullfile(TMP.folder,TMP.name);
    [~,F(i).ete_P] = get_airfoil_surface_pressure_data_from_file(file_name);

    F(i).ic_P = nan*F(i).primal_P;
    TMP = dir(fullfile(folder,folders{i},surface_pat{3}));
    if (isempty(TMP))
        warning('No ic-surface_forces.dat in folder: %s\n',folders{i})
        continue
    end
    file_name = fullfile(TMP.folder,TMP.name);
    [~,F(i).ic_P] = get_airfoil_surface_pressure_data_from_file(file_name);

    F(i).exact_P = nan*F(i).primal_P;
    TMP = dir(fullfile(folder,folders{i},surface_pat{4}));
    if (isempty(TMP))
        warning('No exact-surface_forces.dat in folder: %s\n',folders{i})
        continue
    end
    file_name = fullfile(TMP.folder,TMP.name);
    [~,F(i).exact_P] = get_airfoil_surface_pressure_data_from_file(file_name);
end
end




function H = read_force_histories_from_folder(folder,alpha,reconstruct)
if ( nargin > 2 )&&( reconstruct )
    history_pat = { '*-001-primal_rec-force_history.dat', ...%
                    '*-001-ete_rec-force_history.dat',    ...%
                    '*-001-ic_rec-force_history.dat',     ...%
                    '*-001-exact_rec-force_history.dat' };
else
    history_pat = { '*-001-primal-force_history.dat', ...%
                    '*-001-ete-force_history.dat',    ...%
                    '*-001-ic-force_history.dat',     ...%
                    '*-001-exact-force_history.dat' };
end
H = struct();

% first check for the number of primal-force_history files
% (these will still be output, even if the iterative correction fails) 
TMP = dir(fullfile(folder,'**\',history_pat{1}));
% for now, we will assume a single file per folder
N_files = length(TMP);
if (N_files==0)
    error('no matching primal files in folder')
end

folders = cell(N_files,1);
% get the corresponding folders and primal force data
for i = 1:N_files
    folders{i} = replace(TMP(i).folder,{folder,history_pat{1}},{'',''});
    % fprintf('found folder: %s\n',folders{i});
    H(i).N  = cellfun(@str2double,regexp(folders{i},'\d*','match'));
    H(i).N  = sqrt(prod(H(i).N));
    file_name = fullfile(TMP(i).folder,TMP(i).name);
    H(i).primal_CX = nan(1,1);
    H(i).primal_CY = nan(1,1);
    H(i).primal_CL = nan(1,1);
    H(i).primal_CD = nan(1,1);
    [primal_CX,primal_CY,primal_CL,primal_CD] = ...
        get_force_history_data_from_file(file_name,alpha);
    H(i).primal_CX(:,:) = primal_CX;
    H(i).primal_CY(:,:) = primal_CY;
    H(i).primal_CL(:,:) = primal_CL;
    H(i).primal_CD(:,:) = primal_CD;
end

% get maximum number of successful iterative corrections
max_ic = 1;
for i = 1:N_files
    [ete_CX,~,~,~] = get_force_history_data_from_file(file_name,alpha);
    max_ic = max(max_ic,size(ete_CX,1));
end

% Now loop through the folders and grab the corresponding files
for i = 1:N_files

    H(i).ete_CX = nan(2,1);
    H(i).ete_CY = nan(2,1);
    H(i).ete_CL = nan(2,1);
    H(i).ete_CD = nan(2,1);
    TMP = dir(fullfile(folder,folders{i},history_pat{2}));
    if (isempty(TMP))
        warning('No ete-force_history.dat in folder: %s\n',folders{i})
        continue
    end
    file_name = fullfile(TMP.folder,TMP.name);
    [ete_CX,ete_CY,ete_CL,ete_CD] = ...
        get_force_history_data_from_file(file_name,alpha);
    for j = 1:size(ete_CX,1)
        H(i).ete_CX(j,:) = ete_CX(j);
        H(i).ete_CY(j,:) = ete_CY(j);
        H(i).ete_CL(j,:) = ete_CL(j);
        H(i).ete_CD(j,:) = ete_CD(j);
    end

    H(i).ic_CX = nan(max_ic,1);
    H(i).ic_CY = nan(max_ic,1);
    H(i).ic_CL = nan(1,1);
    H(i).ic_CD = nan(1,1);
    TMP = dir(fullfile(folder,folders{i},history_pat{3}));
    if (isempty(TMP))
        warning('No ic-force_history.dat in folder: %s\n',folders{i})
        continue
    end
    file_name = fullfile(TMP.folder,TMP.name);
    [ic_CX,ic_CY,ic_CL,ic_CD] = ...
        get_force_history_data_from_file(file_name,alpha);
    for j = 1:size(ic_CX,1)
        H(i).ic_CX(j,:) = ic_CX(j);
        H(i).ic_CY(j,:) = ic_CY(j);
        H(i).ic_CL(j,:) = ic_CL(j);
        H(i).ic_CD(j,:) = ic_CD(j);
    end


    H(i).exact_CX = nan(1,1);
    H(i).exact_CY = nan(1,1);
    H(i).exact_CL = nan(1,1);
    H(i).exact_CD = nan(1,1);
    TMP = dir(fullfile(folder,folders{i},history_pat{4}));
    if (isempty(TMP))
        warning('No exact-force_history.dat in folder: %s\n',folders{i})
        continue
    end
    file_name = fullfile(TMP.folder,TMP.name);
    [exact_CX,exact_CY,exact_CL,exact_CD] = ...
        get_force_history_data_from_file(file_name,alpha);
    H(i).exact_CX = exact_CX;
    H(i).exact_CY = exact_CY;
    H(i).exact_CL = exact_CL;
    H(i).exact_CD = exact_CD;
end
end