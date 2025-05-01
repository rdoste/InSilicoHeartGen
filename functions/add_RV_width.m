%     Automated pipeline for large-scale cardiac in silico trials 
%     Copyright (C) 2024 Ruben Doste. Contact: ruben.doste@gmail.com
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>

function add_RV_width(format,width)

    global case_number;
    
%extrude RV
    [ptoRV,facesRV]=read_ply(strcat(format,num2str(case_number),'.ply'));   %surface mesh with no topological defects (can be the same as labels)
    
    TR_RV=triangulation(facesRV,ptoRV);
    normalRV = vertexNormal(TR_RV);
    ptoRV_epi = ptoRV + width * normalRV; %extrude 3 mm
    
%save RV epi
    save_ply(ptoRV_epi,facesRV-1,strcat('RV_epi_',num2str(case_number)));

end

