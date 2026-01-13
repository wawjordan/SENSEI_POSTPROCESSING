function [hfig,DE] = parse_and_plot_new2(dim,r_fac,folders,vars,var_mask,norms,iterations,layers,tags,fmt,color_spec,legend_flag)
if(nargin<11)
    legend_flag = false;
end

N = length(vars);
if isscalar(folders)
    folders = repmat(folders,1,N);
end
if isscalar(var_mask)
    var_mask = repmat(var_mask,1,N);
end
if isscalar(norms)
    norms = repmat(norms,1,N);
end
if isscalar(iterations)
    iterations = repmat(iterations,1,N);
end
if isscalar(layers)
    layers = repmat(layers,1,N);
end
if isscalar(tags)
    tags = repmat(tags,1,N);
end
if isscalar(fmt)
    fmt = repmat(fmt,1,N);
end
if isscalar(color_spec)
    color_spec = repmat(color_spec,1,N);
end
if (length(vars) ~= N || length(var_mask) ~= N || length(iterations) ~= N || length(norms) ~= N || length(layers) ~= N || length(tags) ~= N || length(fmt) ~= N )
    error('inputs of incompatible sizes')
end

% argument cleanup
for i = 1:N
    if ~any( vars(i) == 3:6 )
        layers{i}     = [];
    end
    if ~any( vars(i) == [4,6] )
        iterations{i} = [];
    end
    if size(color_spec{i},1)==1
        color_spec{i} = repmat(color_spec{i},4,1);
    end
end

var_pat = {};

% some old format that I don't remember (keep it for backward compat.)
var_pat{1} = { '*-rec-err_primitive.dat',        ...% reconstruction error     (1)
            '*-te_err.dat',                   ...% truncation error (error) (2)
            '*-de_primal_primitive_1.dat',    ...% Base DE (primitive vars) (3)
            '*-de_corrected_primitive_1.dat', ...% Cor. DE (primitive vars) (4)
            '*-de_primal_conserved_1.dat',    ...% Base DE (conserved vars) (5)
            '*-de_corrected_conserved_1.dat'};   % Cor. DE (conserved vars) (6)
        
% iterative correction + boundary offsets
var_pat{2} = { '*-rec-err_primitive.dat',        ...% reconstruction error     (1)
            '*-te_err.dat',                   ...% truncation error (error) (2)
            '*-de_primal_primitive_offset.dat',    ...% Base DE (primitive vars) (3)
            '*-de_corrected_primitive_ic_offset.dat', ...% Cor. DE (primitive vars) (4)
            '*-de_primal_conserved_offset.dat',    ...% Base DE (conserved vars) (5)
            '*-de_corrected_conserved_ic_offset.dat'};   % Cor. DE (conserved vars) (6)
        
% iterative correction (no offsets)
var_pat{3} = { '*-rec-err_primitive.dat',        ...% reconstruction error     (1)
            '*-te_err.dat',                   ...% truncation error (error) (2)
            '*-de_primal_primitive.dat',    ...% Base DE (primitive vars) (3)
            '*-de_corrected_primitive_ic.dat', ...% Cor. DE (primitive vars) (4)
            '*-de_primal_conserved.dat',    ...% Base DE (conserved vars) (5)
            '*-de_corrected_conserved_ic.dat'};   % Cor. DE (conserved vars) (6)
        
% iterative correction (no offsets)
var_pat{4} = { '*-rec-err_primitive.dat',        ...% reconstruction error     (1)
            '*-te_err.dat',                   ...% truncation error (error) (2)
            '*-de_primal_primitive.dat',    ...% Base DE (primitive vars) (3)
            '*-de_corrected_primitive_ic.dat', ...% Cor. DE (primitive vars) (4)
            '*-de_primal_conserved.dat',    ...% Base DE (conserved vars) (5)
            '*-de_corrected_conserved_ic.dat'};   % Cor. DE (conserved vars) (6)

% boundary offsets (no iterative correction)
var_pat{5} = { '*-rec-err_primitive.dat',        ...% reconstruction error     (1)
            '*-te_err.dat',                   ...% truncation error (error) (2)
            '*-de_primal_primitive_offset.dat',    ...% Base DE (primitive vars) (3)
            '*-de_corrected_primitive_offset.dat', ...% Cor. DE (primitive vars) (4)
            '*-de_primal_conserved_offset.dat',    ...% Base DE (conserved vars) (5)
            '*-de_corrected_conserved_offset.dat'};   % Cor. DE (conserved vars) (6)
        
% boundary offsets (no iterative correction) (compat. without offsets for primal)
var_pat{6} = { '*-rec-err_primitive.dat',        ...% reconstruction error     (1)
            '*-te_err.dat',                   ...% truncation error (error) (2)
            '*-de_primal_primitive.dat',    ...% Base DE (primitive vars) (3)
            '*-de_corrected_primitive_offset.dat', ...% Cor. DE (primitive vars) (4)
            '*-de_primal_conserved.dat',    ...% Base DE (conserved vars) (5)
            '*-de_corrected_conserved_offset.dat'};   % Cor. DE (conserved vars) (6)
        
% no boundary offsets and no iterative correction
var_pat{7} = { '*-rec-err_primitive.dat',        ...% reconstruction error     (1)
            '*-te_err.dat',                   ...% truncation error (error) (2)
            '*-de_primal_primitive.dat',    ...% Base DE (primitive vars) (3)
            '*-de_corrected_primitive.dat', ...% Cor. DE (primitive vars) (4)
            '*-de_primal_conserved.dat',    ...% Base DE (conserved vars) (5)
            '*-de_corrected_conserved.dat'};   % Cor. DE (conserved vars) (6)
        
var_txt = {'rec. error ', 'TE error ','primal DE ','corr. DE ','primal DE ','corr. DE '};
norm_fmt = {'Norm: $L_{1}$','Norm: $L_{2}$','Norm: $L_{\infty}$'};
% base_legend_text_primitive = {'$\\rho$ (%s)','$u$ (%s)','$v$ (%s)','$p$ (%s)'};
% base_legend_text_conserved = {'$\\rho$ (%s)','$\\rho{u}$ (%s)','$\\rho{v}$ (%s)','$\\rho{e}$ (%s)'};
base_legend_text_primitive = {'$\\rho$ %s','$u$ %s','$v$ %s','$p$ %s'};
base_legend_text_conserved = {'$\\rho{\\,\\,\\,}$ %s','$\\rho{u}$ %s','$\\rho{v}$ %s','$\\rho{e}$ %s'};

tmp_mask = sum(cell2mat(var_mask.'),1)>0;
N_plot_var = sum(tmp_mask);

%% string stuff
if legend_flag
    i = 1;
    if (vars(i) == 1)||(vars(i) == 3)||(vars(i) == 4)
        legend_text = {'$\rho$','$u$','$v$','$p$'};
        legend_text = legend_text(tmp_mask);
        if N_plot_var>1
            if N_plot_var<N
                for j = 1:N-N_plot_var
                    legend_text = [legend_text,{''}];
                end
            end
        end
    elseif (vars(i) == 2)||(vars(i) == 5)||(vars(i) == 6)
        legend_text = {'$\rho$','$\rho{u}$','$\rho{v}$','$\rho{e}$'};
    end
    type_text = cell(1,N);
    for i = 1:N
        type_text{i} = [var_txt{vars(i)},tags{i}];
    end
    legend_text = [legend_text,type_text];
else
    type_text = cell(1,N);
    for i = 1:N
        type_text{i} = [var_txt{vars(i)},tags{i}];
    end

    legend_text = cell(4,N);
    for i = 1:N
        if (vars(i) == 1)||(vars(i) == 3)||(vars(i) == 4)
            legend_text(:,i) = cellfun(@(str)sprintf(str,type_text{i}),base_legend_text_primitive,'UniformOutput',false);
        elseif (vars(i) == 2)||(vars(i) == 5)||(vars(i) == 6)
            legend_text(:,i) = cellfun(@(str)sprintf(str,type_text{i}),base_legend_text_conserved,'UniformOutput',false);
        end
    end
end
legend_text = legend_text(:)';
if ~legend_flag
    legend_text = legend_text(logical([var_mask{:}]'));
end

%% read data in
DE = struct();
for i = 1:N
    
    for j = 1:length(var_pat)
        try
            break_loop = true;
            DE(i).DE = DE_parsing(folders{i},vars(i),var_pat{j},dim,r_fac);
        catch
            break_loop = false;
            warning('var_pat{%d} failed',j)
            continue
        end
        if ( break_loop )
            break;
        end
    end
end

%% plotting
lim1 = 10^( floor( log10( DE(1).DE(1).N2   ) ) );
lim2 = 10^( ceil(  log10( DE(1).DE(end).N2 ) ) );

hfig = stdplot(1);
subplot(1,2,1)
hold on
if legend_flag
    p = plot_all_vars_legend_flag(gca,tmp_mask);
    if (sum(tmp_mask)>1)
        if N_plot_var<N
            for j = 1:N-N_plot_var
                ptmp = plot(gca,nan,nan,'LineStyle','none');
                p = [p,{ptmp}];
            end
        end
    end
    ptmp = plot_all_vars_legend_flag_line(gca,fmt);
    p = [p,ptmp];
end
for i = 1:N
    plot_all_vars_DE(gca,DE(i).DE,norms(i),iterations{i},layers(i),var_mask{i},fmt{i},color_spec{i});
end
xlim([lim1,lim2])
% xlim([1,1000])
% ylim([1e-8,1e-3])
% ylim([1e-7,1e-1])

xlabel('$N_{cells}^{1/2}$','Interpreter','latex');
ylabel('Discretization Error Norm','Interpreter','latex')
if (all(norms==norms(1)))
    label_norms(gca,norm_fmt{norms(1)})
end
if legend_flag
    if ( sum(tmp_mask)>1)
        leg1 = legend([p{:}],legend_text,'Interpreter','latex','Location','southwest','NumColumns',2);
    else
        leg1 = legend([p{:}],legend_text,'Interpreter','latex','Location','southwest','NumColumns',1);
    end
    leg1.ItemTokenSize = [18,0];
else
    leg1 = legend(legend_text,'Interpreter','latex','Location','southwest','NumColumns',N);
    leg1.ItemTokenSize = [18,0];
end

subplot(1,2,2)
hold on
if legend_flag
    p = plot_all_vars_legend_flag(gca,tmp_mask);
    if (sum(tmp_mask)>1)
        if N_plot_var<N
            for j = 1:N-N_plot_var
                ptmp = plot(gca,nan,nan,'LineStyle','none');
                p = [p,{ptmp}];
            end
        end
    end
    ptmp = plot_all_vars_legend_flag_line(gca,fmt);
    p = [p,ptmp];
end
for i = 1:N
    plot_all_vars_OOA(gca,DE(i).DE,norms(i),iterations{i},layers(i),var_mask{i},fmt{i},color_spec{i});
end
xlim([lim1,lim2])


% xlim([1,1000])
% ylim([0,3.5])

xlabel('$N_{cells}^{1/2}$','Interpreter','latex');
ylabel('Observed Order of Accuracy','Interpreter','latex')


if legend_flag
    leg2 = legend([p{:}],legend_text,'Interpreter','latex','Location','southwest','NumColumns',2);
    leg2.ItemTokenSize = [22,0];
else
    leg2 = legend(legend_text(),'Interpreter','latex','Location','southwest','NumColumns',N);
    leg2.ItemTokenSize = [22,0];
end

% leg2 = legend(legend_text,'Interpreter','latex','Location','southwest','NumColumns',N);
% leg2 = legend(legend_text,'Interpreter','latex','Location','southwest','FontSize',14,'NumColumns',N);
% leg2.ItemTokenSize = [22,0];
% leg2.ItemTokenSize = [18,0];
if (all(norms==norms(1)))
    label_norms(gca,norm_fmt{norms(1)})
end

end

function DE = DE_parsing(FOLDER,var,var_pat,dim,r_fac)
    DE = read_err_norm_folder(FOLDER, var_pat{var}, dim, r_fac);
    fclose('all');
end

function label_norms(ax,txt)

xpos = 0.05;
ypos = 0.97;

text(ax,xpos,ypos,txt,'Interpreter','latex','Units','Normalized',...
                      'VerticalAlignment','top')

end

function p = plot_all_vars_legend_flag(ax,mask)
hold on
vars = [1,2,3,5];
colors = lines(length(vars));
markers = {'o','s','d','^'};
cnt = 0;
for i = 1:length(vars)
    if ~mask(i)
        continue
    end
    cnt = cnt + 1;
    p{cnt} = plot(ax,nan,nan,'color',colors(i,:));
end
cnt = 0;
for i = 1:length(vars)
    if ~mask(i)
        continue
    end
    cnt = cnt + 1;
    p{cnt}.Marker=markers{i};
    p{cnt}.MarkerEdgeColor=colors(i,:);
    p{cnt}.MarkerSize=3;
    p{cnt}.LineStyle='none';
end
end

function alpha = alpha_calc(x)
min_alpha = 0.1;
max_alpha = 0.25;
alpha = min_alpha + (max_alpha-min_alpha)*sqrt(x);
% alpha = min_alpha + (1-min_alpha)*(2*x-1).^2;
end

function p = plot_all_vars_legend_flag_line(ax,linestyles)
hold on
for i = 1:length(linestyles)
    p{i} = plot(ax,nan,nan,'color','k');
end
for i = 1:length(linestyles)
    p{i}.LineStyle=linestyles{i};
end
end

function p = plot_all_vars_DE(ax,DE,norm,iterations,layers,var_mask,linspec,colorspec)
hold on
vars = [1,2,3,5];
if nargin>6 && size(colorspec,1)==length(vars)
    colors = colorspec;
else
    colors = lines(length(vars));
end
markers = {'o','s','d','^'};
cnt = 0;
for i = 1:length(vars)
    if ~(var_mask(i)), continue; end
    cnt = cnt + 1;
    for j = 1:length(iterations)
        if isscalar(iterations)
            continue;
        else
            alpha = alpha_calc( j/length(iterations) );    
        end
        [ptmp,Ntmp,DEtmp] = plot_DE(ax,DE,vars(i),norm,iterations(j)+1,layers{1},{linspec,'color',[colors(i,:),alpha],'HandleVisibility','off'});
        scatter(Ntmp,DEtmp,3,markers{i},'MarkerEdgeColor',colors(i,:),'MarkerEdgeAlpha',alpha,'LineWidth',ptmp.LineWidth,'HandleVisibility','off');
%         ptmp.Marker=markers{i};
%         ptmp.MarkerEdgeColor=colors(i,:);
%         ptmp.MarkerSize=3;
    end
    if (isempty(iterations))
        p{cnt} = plot_DE(ax,DE,vars(i),norm,[],layers{1},{linspec,'color',colors(i,:)});
    else
        p{cnt} = plot_DE(ax,DE,vars(i),norm,iterations(end)+1,layers{1},{linspec,'color',colors(i,:)});
    end
    p{cnt}.Marker=markers{i};
    p{cnt}.MarkerEdgeColor=colors(i,:);
    p{cnt}.MarkerSize=3;
end
end

function p = plot_all_vars_OOA(ax,DE,norm,iterations,layers,var_mask,linspec,colorspec)
hold on
vars = [1,2,3,5];
if nargin>6 && size(colorspec,1)==length(vars)
    colors = colorspec;
else
    colors = lines(length(vars));
end
markers = {'o','s','d','^'};
cnt = 0;
for i = 1:length(vars)
    if ~(var_mask(i)), continue; end
    cnt = cnt + 1;
    for j = 1:length(iterations)
        if isscalar(iterations)
            continue;
        else
            alpha = alpha_calc( j/length(iterations) );    
        end
        [ptmp,Ntmp,OOAtmp] = plot_OOA(ax,DE,vars(i),norm,iterations(j)+1,layers{1},{linspec,'color',[colors(i,:),alpha],'HandleVisibility','off'});
        scatter(Ntmp,OOAtmp,9,markers{i},'MarkerEdgeColor',colors(i,:),'MarkerEdgeAlpha',alpha,'LineWidth',ptmp.LineWidth,'HandleVisibility','off');
%         ptmp.Marker=markers{i};
%         ptmp.MarkerEdgeColor=colors(i,:);
%         ptmp.MarkerSize=3;
    end
    if (isempty(iterations))
        p{cnt} = plot_OOA(ax,DE,vars(i),norm,[],layers{1},{linspec,'color',colors(i,:)});
    else
        p{cnt} = plot_OOA(ax,DE,vars(i),norm,iterations(end)+1,layers{1},{linspec,'color',colors(i,:)});
    end
    p{cnt}.Marker=markers{i};
    p{cnt}.MarkerEdgeColor=colors(i,:);
    p{cnt}.MarkerSize=3;
end
end



function [p,N_tmp,DE_tmp] = plot_DE(ax,DE,var,norm,iteration,layer,additional_args)
DE_tmp = nan(length(DE),1);
N_tmp  = nan(length(DE),1);

has_iterations = isscalar(DE(1).N_iterations);
has_layers     = isscalar(DE(1).N_layers);

request_iteration = isscalar(iteration);
request_layer     = isscalar(layer);

% handling the case of primal results
if (~has_iterations && request_iteration)
    request_iteration = false;
end

invalid_request = (~has_iterations && request_iteration)||(~has_layers && request_layer);

if invalid_request
    error('The DE structure does not appear to have the requested components');
end
    
if request_layer && request_iteration
    for i = 1:length(DE)
        DE_tmp(i) = DE(i).E(var,norm,iteration,layer);
        N_tmp(i)  = DE(i).N2;
    end
elseif request_layer
    if has_iterations
        for i = 1:length(DE)
            DE_tmp(i) = DE(i).E(var,norm,1,layer);
            N_tmp(i)  = DE(i).N2;
        end
    else
        for i = 1:length(DE)
            DE_tmp(i) = DE(i).E(var,norm,layer);
            N_tmp(i)  = DE(i).N2;
        end
    end
elseif request_iteration
    if has_layers
        for i = 1:length(DE)
            DE_tmp(i) = DE(i).E(var,norm,iteration,1);
            N_tmp(i)  = DE(i).N2;
        end
    else
        for i = 1:length(DE)
            DE_tmp(i) = DE(i).E(var,norm,iteration);
            N_tmp(i)  = DE(i).N2;
        end
    end
else
    for i = 1:length(DE)
        DE_tmp(i) = DE(i).E(var,norm,end,1);
        N_tmp(i)  = DE(i).N2;
    end
end

p = plot(ax,N_tmp,DE_tmp,additional_args{:});
set(ax,'Xscale','log')
set(ax,'Yscale','log')
end

function [p,N_tmp,OOA_tmp] = plot_OOA(ax,DE,var,norm,iteration,layer,additional_args)
OOA_tmp = nan(length(DE),1);
N_tmp  = nan(length(DE),1);

has_iterations = isscalar(DE(1).N_iterations);
has_layers     = isscalar(DE(1).N_layers);

request_iteration = isscalar(iteration);
request_layer     = isscalar(layer);

% handling the case of primal results
if (~has_iterations && request_iteration)
    request_iteration = false;
end

invalid_request = (~has_iterations && request_iteration)||(~has_layers && request_layer);

if invalid_request
    error('The DE structure does not appear to have the requested components');
end
    
if request_layer && request_iteration
    for i = 1:length(DE)
        OOA_tmp(i) = DE(i).OOA(var,norm,iteration,layer);
        N_tmp(i)  = DE(i).N2;
    end
elseif request_layer
    if has_iterations
        for i = 1:length(DE)
            OOA_tmp(i) = DE(i).OOA(var,norm,1,layer);
            N_tmp(i)  = DE(i).N2;
        end
    else
        for i = 1:length(DE)
            OOA_tmp(i) = DE(i).OOA(var,norm,layer);
            N_tmp(i)  = DE(i).N2;
        end
    end
elseif request_iteration
    if has_layers
        for i = 1:length(DE)
            OOA_tmp(i) = DE(i).OOA(var,norm,iteration,1);
            N_tmp(i)  = DE(i).N2;
        end
    else
        for i = 1:length(DE)
            OOA_tmp(i) = DE(i).OOA(var,norm,iteration);
            N_tmp(i)  = DE(i).N2;
        end
    end
else
    for i = 1:length(DE)
        OOA_tmp(i) = DE(i).OOA(var,norm,end,1);
        N_tmp(i)  = DE(i).N2;
    end
end

p = plot(ax,N_tmp,OOA_tmp,additional_args{:});
set(ax,'Xscale','log')
end


function DE = calc_OOA(DE,r_fac)
N_grids = length(DE);
DE(1).OOA = nan*DE(1).E;
for i = 2:N_grids
    den = max(DE(i).E,eps(1));
    Etmp = DE(i-1).E ./ den;
    DE(i).OOA = log(Etmp)./log(r_fac);
end

end

function DE = read_err_norm_folder(folder,pattern,dim,r_fac)
DE = struct();
TMP = dir(fullfile(folder,'**\',pattern));
N_files = length(TMP);
for i = 1:N_files
    file_name = fullfile(TMP(i).folder,TMP(i).name);
    [dat,N,t,N_iterations,N_layers] = read_err_norm_file(file_name);
    DE(i).N1 = N;
    DE(i).N2 = N^(1/dim);
    DE(i).t  = t;
    DE(i).N_iterations = N_iterations;
    DE(i).N_layers = N_layers;
    DE(i).E = dat;
end

[~,ind] = sort([DE.N1]);
DE = DE(ind);
DE = calc_OOA(DE,r_fac);
end

%%%%%%%%%%%%%%

function [dat,N,t,N_iterations,N_layers] = read_err_norm_file(file_name)

N_layers = [];
N_iterations = [];
    
fid = fopen(file_name,'r');

line1 = fgetl(fid);
line2 = fgetl(fid);
line3 = fgetl(fid);


% Check if lines 1-3 contain a single entry
tmp1 = regexp(line1,'\d*','match');
tmp2 = regexp(line2,'\d*','match');
tmp3 = regexp(line3,'\d*','match');

if (length(tmp1) ~= 1)
    % assume no IC, no offsets
    fclose(fid);
    [dat,N,t] = read_err_norm_file_basic(file_name);
    return;
elseif (length(tmp2) ~= 1)
    % assume no IC, but with offets
    fclose(fid);
    [dat,N,t,N_layers] = read_err_norm_file_with_offsets(file_name);
    return
elseif (length(tmp3) ~= 1)
    % assume IC, but no offets
    fclose(fid);
    [dat,N,t,N_iterations] = read_err_norm_ic_file(file_name);
    return
else
    fclose(fid);
    [dat,N,t,N_layers,N_iterations] = read_err_norm_file_with_ic_and_offsets(file_name);
    return
end

end


function [dat,N,t,N_layers,N_iterations] = read_err_norm_file_with_ic_and_offsets(file_name)
% assume IC with offsets
fid = fopen(file_name,'r');

line1 = fgetl(fid);
line2 = fgetl(fid);
fgetl(fid);


% Check if lines 1-3 contain a single entry
tmp1 = regexp(line1,'\d*','match');
tmp2 = regexp(line2,'\d*','match');

% assume that this integer is the number of iterative corrections
N_iterations = str2double(tmp1{1});

% assume that this integer is the number of excluded layers
N_layers = str2double(tmp2{1});

% the next line should be the (grid size, 1D grid size, time)
line = fgetl(fid);
tmp = sscanf(line,'%f');
N = tmp(1);

has_t = false;
if (length(tmp)==3)
    t = zeros(N_iterations+1,1);
    t(1) = tmp(3);
    has_t = true;
else
    t = [];
end

% this should be the first line with error norms
line = fgetl(fid);
tmp = sscanf(line,'%f');
N_eqns = length(tmp);

% allocate
dat = zeros(N_eqns,3,N_iterations+1,N_layers+1);

% rewind
frewind(fid);
% read the first two lines again
fgetl(fid);
fgetl(fid);
% now loop
for j = 0:N_iterations
    line = fgetl(fid); % (iteration) not needed
    if (ischar(line))
        for k = 0:N_layers
            line = fgetl(fid); % (grid size, 1D grid size, time)
            if has_t
                if (k==0)
                    tmp = sscanf(line,'%f');
                    t(j+1) = tmp(3);
                end
            end
            for i = 1:3 % norms (L1, L2, L_inf)
                line = fgetl(fid);
                tmp = sscanf(line,'%f');
                dat(:,i,j+1,k+1) = tmp;
            end
        end
    else
        break
    end
end

fclose(fid);
end


function [dat,N,t,N_iterations] = read_err_norm_ic_file(file_name)
% assume IC, no offsets
fid = fopen(file_name,'r');

line1 = fgetl(fid);
fgetl(fid);


% Check if lines 1-2 contain a single entry
tmp1 = regexp(line1,'\d*','match');

% assume that this integer is the number of iterative corrections
N_iterations = str2double(tmp1{1});

% the next line should be the (grid size, 1D grid size, time)
line = fgetl(fid);
tmp = sscanf(line,'%f');
N = tmp(1);

has_t = false;
if (length(tmp)==3)
    t = zeros(N_iterations+1,1);
    t(1) = tmp(3);
    has_t = true;
else
    t = [];
end

% this should be the first line with error norms
line = fgetl(fid);
tmp = sscanf(line,'%f');
N_eqns = length(tmp);

% allocate
dat = zeros(N_eqns,3,N_iterations+1);

% rewind
frewind(fid);
% read the first line again
fgetl(fid);
% now loop
for j = 0:N_iterations
    fgetl(fid); % (iteration) not needed
    line = fgetl(fid); % (grid size, 1D grid size, time)
    if has_t
        tmp = sscanf(line,'%f');
        t(j+1) = tmp(3);
    end
    for i = 1:3 % norms (L1, L2, L_inf)
        line = fgetl(fid);
        tmp = sscanf(line,'%f');
        dat(:,i,j+1) = tmp;
    end
end
fclose(fid);
end



function [dat,N,t,N_layers] = read_err_norm_file_with_offsets(file_name)
% assume no IC, offsets
fid = fopen(file_name,'r');

line = fgetl(fid);

tmp = regexp(line,'\d*','match');

% assume that this integer is the number of excluded layers (+1)
N_layers = str2double(tmp{1});

% the next line should be the (grid size, 1D grid size, time)
line = fgetl(fid);
tmp = sscanf(line,'%f');
N = tmp(1);

if (length(tmp)==3)
    t = tmp(3);
else
    t = [];
end

line = fgetl(fid);
tmp = sscanf(line,'%f');
N_eqns = length(tmp);

% allocate
dat = zeros(N_eqns,3,N_layers+1);

% rewind
frewind(fid);
% read the first line again
fgetl(fid);
% now loop
for j = 0:N_layers
    fgetl(fid); % not needed anymore
    for i = 1:3 % remaining norms (L2, L_inf)
        line = fgetl(fid);
        tmp = sscanf(line,'%f');
        dat(:,i,j+1) = tmp;
    end
    
end

fclose(fid);
end

function [dat,N,t] = read_err_norm_file_basic(file_name)
fid = fopen(file_name,'r');

% the next line should be the (grid size, 1D grid size, time)
line = fgetl(fid);
tmp = sscanf(line,'%f');
N = tmp(1);

if (length(tmp)==3)
    t = tmp(3);
else
    t = [];
end

line = fgetl(fid);
tmp = sscanf(line,'%f');
N_eqns = length(tmp);
dat = zeros(N_eqns,3);
dat(:,1) = tmp;

for i = 2:3
    line = fgetl(fid);
    tmp = sscanf(line,'%f');
    dat(:,i) = tmp;
end

fclose(fid);
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