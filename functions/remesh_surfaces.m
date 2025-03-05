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
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.

function  geo4=remesh_surfaces(edge_length)
    
    global case_number; 
    
    
  
    disp('Surf mesh generation');
    
    %Remesh faces
        surfname=["Closed_mesh_epi","Closed_mesh_endo_LV","Closed_mesh_endo_RV"];
        for ind=1:3
            [pto0,face0]=read_ply(strcat(surfname(ind),'.ply'));
            opt.gridsize=1/4;
            opt.elemsize=edge_length;
            [pto2,faces2]=remeshsurf(pto0,face0,opt);
            faces2=faces2(:,1:end-1);
            [pto2,faces2]=meshcheckrepair(pto2,faces2,1);
            [pto2,faces2]=surfreorient(pto2,faces2);
            [pto2,faces2]=meshresample(pto2,faces2,0.98); 
        
            save_ply(pto2,faces2-1,strcat('Remesh_',surfname(ind)));
        end
    
    
    %merge closed geometries
        
         [geo1]=vtkRead('Remesh_Closed_mesh_epi.ply');
         [geo2]=vtkRead('Remesh_Closed_mesh_endo_LV.ply');
         [geo3]=vtkRead('Remesh_Closed_mesh_endo_RV.ply');
        
         geo4 = vtkAppendFilter({geo1,geo2,geo3});
         vtkWrite(geo4,strcat('Closed_final_',num2str(case_number),'.ply'));



end