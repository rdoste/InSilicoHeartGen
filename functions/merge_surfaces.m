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

function labels0=merge_surfaces(Prefix_LV,Prefix_LV_epi,Prefix_RV)

    global case_number

    [node1,faces1]=read_ply(strcat(Prefix_RV,num2str(case_number),'.ply'));   %surface mesh with no topological defects (can be the same as labels)
    [node2,faces2]=read_ply(strcat('RV_epi_',num2str(case_number),'.ply'));   %surface mesh with no topological defects (can be the same as labels)
    [node3,faces3]=read_ply(strcat(Prefix_LV,num2str(case_number),'_mod.ply'));   %surface mesh with no topological defects (can be the same as labels)
    [node4,faces4]=read_ply(strcat(Prefix_LV_epi,num2str(case_number),'.ply'));   %surface mesh with no topological defects (can be the same as labels)
    
    [newnode_epi,newelem_epi]=surfboolean(node2,faces2,'union', node4,faces4);
    save_ply(newnode_epi,newelem_epi-1,'mergedepi');
    [newnode_endo,newelem_endo]=surfboolean(node1,faces1,'union', node3,faces3);
    %save_ply(newnode_endo,newelem_endo-1,'mergedendo');
    
    [newnode_final,newelem_final]=surfboolean(newnode_epi,newelem_epi,'all',newnode_endo,newelem_endo);
    %save_ply(newnode_final,newelem_final-1,'merged_final');
    
    save_ply(node1,faces1-1,'Closed_mesh_endo_RV');
    
    %% initial labelling
        centroid_final=meshcentroid(newnode_final,newelem_final);
    %epi_LV
        centroid_epi_RV=meshcentroid(node2,faces2);
        TR_epi_RV=delaunayTriangulation(centroid_epi_RV);
        [NN_epi_RV,dist_epi_RV]=nearestNeighbor(TR_epi_RV,centroid_final);
        
        centroid_endo_RV=meshcentroid(node1,faces1);
        TR_endo_RV=delaunayTriangulation(centroid_endo_RV);
        [NN_endo_RV,dist_endo_RV]=nearestNeighbor(TR_endo_RV,centroid_final);
        
        centroid_epi_LV=meshcentroid(node4,faces4);
        TR_epi_LV=delaunayTriangulation(centroid_epi_LV);
        [NN_epi_LV,dist_epi_LV]=nearestNeighbor(TR_epi_LV,centroid_final);
        
        centroid_endo_LV=meshcentroid(node3,faces3);
        TR_endo_LV=delaunayTriangulation(centroid_endo_LV);
        [NN_endo_LV,dist_endo_LV]=nearestNeighbor(TR_endo_LV,centroid_final);
        
        label=zeros(1,length(centroid_final));
        [~,label2]=min([dist_endo_LV,dist_epi_LV,dist_endo_RV,dist_epi_RV],[],2);
        
    %write results
        write_vtk_surf('labels0.vtk',newnode_final,newelem_final,label2); %creation of the first mesh with surface labels
        labels0=struct();
        labels0.points=newnode_final;
        labels0.cells=newelem_final;
        labels0.cellTypes=uint8(ones(length(labels0.cells),1)*5);
        labels0.cellData=label2;
end
