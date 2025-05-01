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

function CSVFilesgeneration_f_UKBB(name,matfilename)
%generation of csv files required for running eikonal

%execute inside folder
    
    
    load(matfilename);
    mkdir CSVFiles
    cd("CSVFiles")
    
    %node coordinates
    writematrix(v,strcat(name,'_xyz.csv'));
    
    %Fibers on nodes
    writematrix(F,strcat(name,'_nodefield_fibre.csv'));
    writematrix(F_S,strcat(name,'_nodefield_sheet.csv'));
    writematrix(F_N,strcat(name,'_nodefield_normal.csv'));
    
    %cobi coordinates
    writematrix(a2b_cobi,strcat(name,'_nodefield_cobiveco-ab.csv'));
    writematrix(lvrv_cobi,strcat(name,'_nodefield_cobiveco-tv.csv'));
    writematrix(Tphi3,strcat(name,'_nodefield_cobiveco-tm.csv'));
    writematrix(r,strcat(name,'_nodefield_cobiveco-rt.csv'));
    writematrix([a2b_cobi,r,tm_cobi,lvrv_cobi],strcat(name,'_cobiveco.csv'));
    
    %gradients
    writematrix(a2p,strcat(name,'_nodefield_cobiveco-aprt.csv'));
    writematrix(r2l_geo,strcat(name,'_nodefield_cobiveco-rvlv.csv'));
    writematrix(r2l,strcat(name,'_nodefield_cobiveco-rvlv-projection.csv'));

    %materials tetra
    writematrix(Plug_tetra,strcat(name,'_material_tetra.csv'));
    
    %tetra connectivity
    writematrix(f(:,2:end)+1,strcat(name,'_tetra.csv'));
    
    %boundaries
    face=volface(f(:,2:5)+1);
        %RV vs LV in boundaries
        RV_ind=intersect(find(and(lvrv_cobi==1,Ventricle==1)),unique(face));
        LV_ind=intersect(find(and(lvrv_cobi==0,Ventricle==-2)),unique(face));
    
    writematrix(RV_ind,strcat(name,'_boundarynodefield_ep-rvnodes.csv'));
    writematrix(LV_ind,strcat(name,'_boundarynodefield_ep-lvnodes.csv'));

    %write cobiveco originals
    % writematrix(a2b_projected,strcat(name,'_nodefield_cobiveco-ab-original.csv'));
    % writematrix(tv_projected,strcat(name,'_nodefield_cobiveco-tv-original.csv'));
    % writematrix(tm_projected,strcat(name,'_nodefield_cobiveco-tm-original.csv'));
    % writematrix(rt_projected,strcat(name,'_nodefield_cobiveco-rt-original.csv'));
    % writematrix([a2b_projected,rt_projected,tm_projected,tv_projected],strcat(name,'_cobiveco-original.csv'));
    cd ..

end
