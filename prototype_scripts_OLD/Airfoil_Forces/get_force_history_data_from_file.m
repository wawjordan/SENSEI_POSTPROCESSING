function [CX,CY,CL,CD] = get_force_history_data_from_file(file_name,alpha)
zones = tec2mat_structured_condensed(file_name);
CX = getfield(zones.ZONE,'CX');
CY = getfield(zones.ZONE,'CY');
CL = CY*cosd(alpha) - CX*sind(alpha);
CD = CY*sind(alpha) + CX*cosd(alpha);
end