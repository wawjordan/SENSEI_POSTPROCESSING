function DATA = get_structured_from_sensei_tec2mat(ZONES,zone_num)

ZONE = ZONES(zone_num);
nvars = length(ZONE.VARLOC);
DATA = struct();
cell_var_cnt = 0;
node_var_cnt = 0;

% imax = ZONE.I;
if isfield(ZONE,'J')
    jmax = ZONE.J;
else
    jmax = 1;
end

if isfield(ZONE,'K')
    kmax = ZONE.K;
else
    kmax = 1;
end

o = 0;
for n = 1:nvars
    DATA(n+o).name = ZONE.ZONEVARlist{n};
    if ZONE.VARLOC(n) == 0 % nodal
        node_var_cnt = node_var_cnt + 1;
        if jmax > 1
            if kmax > 1
                DATA(n+o).data = ZONE.data_nodal(:,:,:,node_var_cnt);
            else
                DATA(n+o).data = ZONE.data_nodal(:,:,node_var_cnt);
            end
        else
            DATA(n+o).data = ZONE.data_nodal(:,node_var_cnt);
        end
        % add cell-centered coordinates
        if any( DATA(n+o).name==['X','Y','Z'])
            hi = size(DATA(n+o).data,[1,2,3]);
            o = o + 1;
            DATA(n+o).name = [DATA(n+o-1).name,'C'];
            DATA(n+o).data = zeros(hi(1)-1,hi(2)-1,hi(3)-1);
            for k = 1:hi(3)-1
                for j = 1:hi(2)-1
                    for i = 1:hi(1)-1
                        DATA(n+o).data(i,j,k) = mean(DATA(n+o-1).data(i:i+1,j:j+1,k:k+1),"all");
                    end
                end
            end
        end         
    elseif ZONE.VARLOC(n) == 1 % cell-centered
        cell_var_cnt = cell_var_cnt + 1;
        if jmax > 2
            if kmax > 2
                DATA(n+o).data = ZONE.data_cell(:,:,:,cell_var_cnt);
            else
                DATA(n+o).data = ZONE.data_cell(:,:,cell_var_cnt);
            end
        else
            DATA(n+o).data = ZONE.data_cell(:,cell_var_cnt);
        end
    else
        error('VARLOC should either be 0 or 1')
    end
end

end