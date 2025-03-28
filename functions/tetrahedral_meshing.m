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

function MeshCoarse=tetrahedral_meshing(surf0,elemsize,pt_LV, pt_RV)

     %MeshCoarse--> Resulting volumetric mesh
     %labels_orig_name-->Name of the initial label mesh. Used to find inner
     %points of the RV and LV cavities.
 
     %surf0 --> Surface triangular mesh used as input
     %elemsize-->tetrahedral element length
    
    %% Remesh Tetra
    disp('New mesh generation');
      
    opt.elemsize=elemsize; 
    Elementvol= (opt.elemsize^3/(6*sqrt(2)));
      

    % Surface Remeshing
     [pto1,faces1]=meshcheckrepair(surf0.points,surf0.cells,1);
     [pto1,faces1]=surfreorient(pto1,faces1);
     [pto1,faces1]=meshresample(pto1,faces1,1.0); 
    

    % Volume meshing
    
        [pt,p0,v0,t,idx]=surfinterior(pto1,faces1);
        [node,elem,~]=surf2mesh(pto1,faces1,[-100,0,0],[200,300,200],1.0,Elementvol,pt,[pt_RV;pt_LV]);
        opt.reratio=1.41;       
        elem=elem(:,1:end-1);        
        [node,elem,faces]=meshrefine(node,elem,[],opt);
        [elem, evol]=meshreorient(node,elem);
        MeshCoarse.points=node;
        MeshCoarse.cells=int32(elem);
        MeshCoarse.cellTypes=uint8(ones(length(MeshCoarse.cells),1)*10);
        

end
