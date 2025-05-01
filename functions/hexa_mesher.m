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

function [ALG,tet_ID,bar_final,X,Y,Z]=hexa_mesher(mesh_name,cube_length,centered)

    %creation of hexahedral mesh for monoAlg3D simulations 
    
    %cube_length-->hexahedral edge length
    %centered-->Flag. If centered==1, translate the mesh to center the minimum in coord
    %(0,0,0), (specific for Monoalg3D)
    
    %hexa_mesher('Coarse.vtu',0.04,1);
    
    %ALG--> .alg mesh format (cube center, half_edge_x,half_edge_y,half_edge_z)
    %tet--> Tetrahedra ID where the cube center is located
    %bar--> Barycentric coordinates of the point
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    %caculate limits of the tet mesh
    
    %read original mesh
    MeshCoarse=vtkRead(mesh_name);
    
    if centered==1
        MeshCoarse.points=[MeshCoarse.points(:,1)-min(MeshCoarse.points(:,1)),MeshCoarse.points(:,2)-min(MeshCoarse.points(:,2)),MeshCoarse.points(:,3)-min(MeshCoarse.points(:,3))];
    end
    
    %specify edge_length
    cube_length_half=cube_length/2;
    
    
    %create grid
    [X,Y,Z] = ndgrid(min(MeshCoarse.points(:,1))+cube_length_half:cube_length:max(MeshCoarse.points(:,1)),min(MeshCoarse.points(:,2))+cube_length_half:cube_length:max(MeshCoarse.points(:,2)),min(MeshCoarse.points(:,3))+cube_length_half:cube_length:max(MeshCoarse.points(:,3)));
    
    coords = reshape([X(:), Y(:), Z(:)], [], 3);
    
    
    
    %find if center of cubes are inside tet mesh
    T=triangulation(double(MeshCoarse.cells),double(MeshCoarse.points));
    [tet,bar]=pointLocation(T,double(coords));
    index_inside=~isnan(tet);
    
    tet_ID=tet(index_inside);
    bar_final=bar(index_inside,:);
    
    %transform center cubs into coordinates
     Coordinates_total=[coords(index_inside,1)-cube_length_half,coords(index_inside,2)-cube_length_half,coords(index_inside,3)-cube_length_half,coords(index_inside,1)+cube_length_half,coords(index_inside,2)-cube_length_half,coords(index_inside,3)-cube_length_half,coords(index_inside,1)+cube_length_half,coords(index_inside,2)+cube_length_half,coords(index_inside,3)-cube_length_half,...
                            coords(index_inside,1)-cube_length_half,coords(index_inside,2)+cube_length_half,coords(index_inside,3)-cube_length_half,coords(index_inside,1)-cube_length_half,coords(index_inside,2)-cube_length_half,coords(index_inside,3)+cube_length_half,coords(index_inside,1)+cube_length_half,coords(index_inside,2)-cube_length_half,coords(index_inside,3)+cube_length_half,...                     
                            coords(index_inside,1)+cube_length_half,coords(index_inside,2)+cube_length_half,coords(index_inside,3)+cube_length_half,coords(index_inside,1)-cube_length_half,coords(index_inside,2)+cube_length_half,coords(index_inside,3)+cube_length_half];
     
     Coordinates_total_reshape=single(reshape(Coordinates_total',[3,size(Coordinates_total,1)*size(Coordinates_total,2)/3]))';
     Coordinates=unique((Coordinates_total_reshape),'stable','rows');
        
     [~,Connectivity_pre]=ismember(Coordinates_total_reshape,Coordinates,'rows');
        
     Connectivity=reshape(Connectivity_pre',[8,length(Coordinates_total)])';
    
     %% write results
     %ALG
    
    ALG=single([coords(index_inside,1),coords(index_inside,2),coords(index_inside,3),repmat(cube_length_half,length(coords(index_inside,1)),3)]*1000);

     %VTK
        namehex=mesh_name(1:end-4);
       
        name=strcat('hex_',namehex,'.vtk');
        f=Connectivity-1;
        fileID = fopen(name,'w');
        fprintf(fileID,'# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\n' );
        fprintf(fileID,'POINTS %d',length(Coordinates));
        fprintf(fileID,' float\n');
        fprintf(fileID,'%g %g %g %g %g %g %g %g %g \n',Coordinates');
        fprintf(fileID,'\n');
        fprintf(fileID,'\n');
        fprintf(fileID,'CELLS %d %d\n',length(Connectivity),length(Connectivity)*9);
        fprintf(fileID,'8 %d %d %d %d %d %d %d %d\n',f');
        fprintf(fileID,'\n');
        fprintf(fileID,'CELL_TYPES %d\n',length(Connectivity));
        fprintf(fileID,'%d \n',zeros(size(f(:,1)))+12);
        fprintf(fileID,'\n');
        fprintf(fileID,'POINT_DATA %d\n',length(Coordinates));
        fclose(fileID);  

     toc 

end