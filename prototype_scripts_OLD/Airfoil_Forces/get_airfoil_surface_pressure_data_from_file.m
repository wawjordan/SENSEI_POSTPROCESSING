function [XC,P] = get_airfoil_surface_pressure_data_from_file(file_name)
zones = tec2mat_structured_condensed(file_name);
Nzones = numel(zones);
XC = getfield(zones(1).ZONE,'XC');
varname = 'CP';
if (isfield(zones(1).ZONE,'P'))
    varname = 'P';
end
[Ni,Nj] = size(getfield(zones(1).ZONE,varname));
P = zeros(Ni,Nj,Nzones);
for n = 1:Nzones
    P(:,:,n) = getfield(zones(n).ZONE,varname);
end
P = squeeze(P);
end