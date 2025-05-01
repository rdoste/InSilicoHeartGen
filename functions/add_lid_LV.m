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

function add_lid_LV()
        [~,pto0,~,car0,label]=read_labels2('labels0.vtk'); %mesh with labels
        car=car0+1;
        
        TR_lab=triangulation(double(car),pto0);
        V = vertexAttachments(TR_lab,double(unique(car)));
        index_labelp=cellfun(@(v)v(1),V);
        labelp=zeros(length(pto0),1);
        labelp(unique(car))=label(index_labelp);
        % write_vtk_surf('labelsp.vtk',pto0,car,labelp);
        
        
        newelem_final=car;
        newnode_final=pto0;
    
    %locate boundary nodes
        newelem_final1 = newelem_final(label==1,:);
        newelem_final2 = newelem_final(label~=1,:);
        boundary1=compute_boundary(newelem_final1);
        boundary2=compute_boundary(newelem_final2);
        b1=newnode_final(boundary1,:);
        b2=newnode_final(boundary2,:);
        
         Chull1 = convhull (b1);
         Chull2 = convhull (b2);

    
    
    %merge with mesh
    
        newelem_final0=[newelem_final(label==2,:);newelem_final(label==4,:)];
        [newnode_closed,newelem_closed]=surfboolean(b2,Chull2,'all', newnode_final,newelem_final0);
        save_ply(newnode_closed,newelem_closed-1,'Closed_mesh_epi');
        
        
        [newnode_closed,newelem_closed]=surfboolean(b1,Chull1,'all', newnode_final,newelem_final1);
        save_ply(newnode_closed,newelem_closed-1,'Closed_mesh_endo_LV');
        
     

end