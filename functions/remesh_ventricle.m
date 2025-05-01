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

function [newmesh_points,newmesh_faces]=remesh_ventricle(points,faces,elementsize,downsamplingrate,scale,meshtype,writemesh)
%tag 'open' or 'closed' biventricular (open without valves, closed with
%valve plane)
%Create  tetrahedral mesh



% check scale size
% is in cm, needs to be changed to mm
% scale=1;



%%

%Remesh Tetra
 facecell=finddisconnsurf(faces); %extract epi, endoLV and endoRV
%subdivide in elements

disp('Surf mesh generation');

%Remesh faces
if meshtype=="closed"
    surfname=["surf1","surf2","surf3"];
    end_ind=3;
elseif meshtype=="UKBB"
    surfname=["Remesh_Closed_mesh_epi_LV","Remesh_Closed_mesh_endo_LV","Remesh_Closed_mesh_endo_RV"];
    end_ind=3;
elseif meshtype=="open" || meshtype=="cut" || meshtype=="cut_generic"
    end_ind=1;
    surfname=["Remeshed_ventricle",0,0];
else
    error('wrong tag')
end

    for ind=1:end_ind
    
        opt.gridsize=1/4;
        opt.elemsize=elementsize;
        face2=facecell(ind);
        [points2,faces2]=remeshsurf(points.*scale,face2{1,1},opt);
        faces2=faces2(:,1:end-1);
        [points2,faces2]=meshcheckrepair(points2,faces2,1);
        [points2,faces2]=surfreorient(points2,faces2);
        [points2,faces2]=meshresample(points2,faces2,downsamplingrate);
        %save_ply(pto2,faces2-1,surfname(ind));
        coarsemesh_points0{ind}=points2./scale;
        coarsemesh_faces0{ind}=faces2;
    end

    if end_ind==2
     [newmesh_points,newmesh_faces]=surfboolean(coarsemesh_points0{1},coarsemesh_faces0{1},'all',coarsemesh_points0{2},coarsemesh_faces0{2});
    elseif end_ind==3
        [newmesh_points,newmesh_faces]=surfboolean(coarsemesh_points0{1},coarsemesh_faces0{1},'all',coarsemesh_points0{2},coarsemesh_faces0{2},'all',coarsemesh_points0{3},coarsemesh_faces0{3});
    elseif end_ind==1
        newmesh_points=coarsemesh_points0{1};
        newmesh_faces=coarsemesh_faces0{1};

    else
        error('mesh boolean operation that involves more than 3 meshes needs to be added in the code')
    end

    %write
    if writemesh==1
        geo=struct();
        geo.points=newmesh_points;
        geo.cells=int32(newmesh_faces);
        geo.cellTypes=uint8(ones(length(geo.cells),1)*5);
        vtkWrite(geo,strcat('Remeshed.ply'));

    end
end