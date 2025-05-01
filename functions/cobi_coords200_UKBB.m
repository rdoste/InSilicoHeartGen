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
function [r, a2b_cobi_coord_final,tm_cobi,lvrv_cobi,apex_cobi_id]=cobi_coords200_UKBB(v,labelp,a2b,a2p,T_sept_cobi,elem,label,labelf3,face,face3,Fid,Ventricle,Tv_cobi,L,GT)
%calculation of the cobiveco coordinates

index_interior=setdiff(1:length(Ventricle),(unique(face)))';
labelp(index_interior)=0;
%Ridge ant vs post according apex a2p coordinate

if max(Fid)==5 || max(Fid)==6
 apex_nodes=find(labelp==5);
elseif max(Fid)>=18 
 apex_nodes=find(labelp==12);
end
 apex_node=apex_nodes(1);

 C_apex=a2p(apex_node);


 %% Obtain surfaces that intersect between septum and rest of the heart
    T_sept_cobi(abs(Ventricle)==3)=0;
    body=T_sept_cobi(elem);

 if max(Fid==18)
    body=T_sept_cobi(elem);
    bodya2b=min(a2b(elem),[],2);

    %removal of the heart base
    %locate max point using a2b gradient
    %minimum valve point in a2b coordinate
    %valve points
    indValve=(labelp==9 | labelp==10) |(labelp==13 | labelp==14);
    %indValve=(labelp==9);

    Threshold_val=min(a2b(indValve));
    Threshold_val=Threshold_val-0.05;
    %elements above Threshold
    index_heart_elem=bodya2b<(Threshold_val);

    
    %faces intersection
    body_bound_min=find(min(body,[],2)<0.5 & bodya2b<Threshold_val);
    body_bound_max=find(max(body,[],2)>=0.5 & bodya2b<Threshold_val);
    body_bound_intersect=intersect(body_bound_min,body_bound_max);

else
              
    %faces with elements 0 and 360
    body_bound_min=find(min(body,[],2)<0.5);
    body_bound_max=find(max(body,[],2)>=0.5);
    body_bound_intersect=intersect(body_bound_min,body_bound_max);

    index_heart_elem=logical(ones(size(body(:,1))));


end
    
    %calculate volfaces
    surf_intersect=volface(elem(body_bound_intersect,:)); 
    surf_min=volface(elem(body_bound_min,:)); 
    surf_max=volface(elem(body_bound_max,:)); 
    
    %bondary surfaces
    surf_1=intersect(surf_intersect,surf_max,'rows');
    surf_0=intersect(surf_intersect,surf_min,'rows');
    
    repeated=intersect(surf_1,surf_0,'rows');
    
    surf_1=setdiff(surf_1,repeated,'rows');
    surf_0=setdiff(surf_0,repeated,'rows');

    %Anterior vs posterior Ridge (I take the surf_0 as reference)

    surf_Ridge_a2p=min(a2p(surf_0),[],2);
    surf_a2p_tags=ones(size(surf_Ridge_a2p))*40;
    surf_a2p_tags(surf_Ridge_a2p<=C_apex)=41;

 
    % Repeat for the outer surface used in the Vfree boundary conditions
    surf_Ridge_a2p_v2=min(a2p(surf_1),[],2);
    surf_a2p_tags_v2=ones(size(surf_Ridge_a2p_v2))*42;
    surf_a2p_tags_v2(surf_Ridge_a2p_v2<=C_apex)=43;

   %force that there are only 2 clusters 

   for ind=40:41
     cluster=finddisconnsurf(surf_0(surf_a2p_tags==ind,:));
       while size(cluster,2)>1
        [~, max_index] = max(cellfun('size', cluster, 1));%detect maximum size array
        cluster{:,max_index}=[];%remove the biggest surface
        labelsf_isolated =ismember(surf_0,cat(1,   cluster{:}),'rows');%isolated faces
        vals=surf_a2p_tags(labelsf_isolated);
            if vals(1)==40
                surf_a2p_tags(labelsf_isolated)=41;
            else
                surf_a2p_tags(labelsf_isolated)=40;
            end
        cluster=finddisconnsurf(surf_0(surf_a2p_tags==ind,:));
       end
   end

   for ind=42:43
     cluster=finddisconnsurf(surf_1(surf_a2p_tags_v2==ind,:));
       while size(cluster,2)>1
        [~, max_index] = max(cellfun('size', cluster, 1));%detect maximum size array
        cluster{:,max_index}=[];%remove the biggest surface
        labelsf_isolated =ismember(surf_1,cat(1,   cluster{:}),'rows');%isolated faces
        vals=surf_a2p_tags_v2(labelsf_isolated);
            if vals(1)==42
                surf_a2p_tags_v2(labelsf_isolated)=43;
            else
                surf_a2p_tags_v2(labelsf_isolated)=42;
            end
        cluster=finddisconnsurf(surf_1(surf_a2p_tags_v2==ind,:));
       end
   end


    %write new labels
    labelR=[label;surf_a2p_tags];               
    faceR= [face;surf_0];
    %write_vtk_surf('labelsR.vtk',v,faceR,labelR);
    
    labelR2=[label;surf_a2p_tags_v2];     
    faceR2=[face;surf_1];
    %write_vtk_surf('labelsR2.vtk',v,faceR2,labelR2);

    %distance calculations septum
    
    idx_sept_elems=(ismember(1:length(elem),setdiff(body_bound_max,body_bound_intersect)));
    Sept_elem=sortconnectivities(elem(and(idx_sept_elems',index_heart_elem),:));
    Sept_points=unique(elem(and(idx_sept_elems',index_heart_elem),:));


    GSept = grad(double(v(Sept_points,:)), double(Sept_elem));
    %GFree = grad(double(v(Free_points,:)), double(Free_elem));

    %distance calculations free vol
    % idx_sept_elems2=(ismember(1:length(elem),body_bound_max));
    Free_elem=sortconnectivities(elem(and(~idx_sept_elems',index_heart_elem),:));
    Free_points=unique(elem(and(~idx_sept_elems',index_heart_elem),:));
    GFree = grad(double(v(Free_points,:)), double(Free_elem));



    %detect points in Ridge
  
    indx_post_row=(labelR==41);
    indx_post_row=unique(faceR(indx_post_row,:));

    indx_ant_row=(labelR==40);
    indx_ant_points=unique(faceR(indx_ant_row,:));
    C_apex_points_id=intersect( indx_ant_points,indx_post_row);
    indx_ant_points=setdiff( indx_ant_points,indx_post_row);

    indx_post_row_free=(labelR==41);
    indx_post_row_free=unique(faceR(indx_post_row_free,:));

    indx_ant_row_free=(labelR==40);
    indx_ant_points_free=unique(faceR(indx_ant_row_free,:));
    indx_ant_points_free=setdiff( indx_ant_points_free,indx_post_row_free);
    
%% Prepare tm, r and a2b coord
% L=cotmatrix(v, elem);
% GT= grad(double(v), double(elem));

%% Transventricular coordiante
lvrv_cobi=ones(size(Ventricle));
lvrv_cobi(Tv_cobi<=0.5)=0;

%% Transmural coordinate

if max(Fid)==5 || max(Fid)==6
    Ventricle(labelp==4)=4;
    idsEpi=find(abs(Ventricle)==3 );
    idsEndo=find(Ventricle==1 | Ventricle==-2);
    elseif max(Fid)==18
    idsEpi=find(abs(Ventricle)==3 );
    idsEndo=find(Ventricle==1 | Ventricle==-2);
end

%septal points
idsSept_f=(labelf3==30); %septal faces
idsSept=unique(face3(idsSept_f,:));
idsSept=setdiff(idsSept,idsEpi);
idsSept=setdiff(idsSept,idsEndo);

ids = [idsEpi; idsSept; idsEndo];
val = [zeros(size(idsEpi,1)+size(idsSept,1),1); ones(size(idsEndo))];
tmLaplace = solveLaplace(L, ids, val, 1e-8, 5000);

T = normalizedGradField(GT, tmLaplace,1e-9, false, v, elem);

tmDistEpi = solveTrajectDist(GT, T, [idsEpi; idsSept], zeros(size(idsEpi,1)+size(idsSept,1),1), 1e-8, 3000);
tmDistEndo = solveTrajectDist(GT, -T, idsEndo, zeros(size(idsEndo)), 1e-8, 3000);

tm_cobi = tmDistEpi./(tmDistEpi+tmDistEndo);
tm_cobi=min(max(tm_cobi,0),1);

tmFlipped = tm_cobi;
tmFlipped(lvrv_cobi==1) = -tmFlipped(lvrv_cobi==1);
tmGrad_c = normalizedGradField(GT, tmFlipped,  1e-8, false, v, elem);


%% a2b direction   
% to do  (calculate until the threshold)
% cobi apex for the a2b direction

TR3=triangulation(double(elem),double(v));
indx_C_apex_relative=(ismember(Sept_points,C_apex_points_id));
tm_cobi_Sept=tm_cobi(Sept_points);
[~,idx_max]=min(tm_cobi_Sept(indx_C_apex_relative));%get the most centric of the C_apex points
center_apex=v(C_apex_points_id(idx_max),:);

[apex_cobi_id,~]=nearestNeighbor(TR3,center_apex); %ID of the apex point located in the apex
% basal points id
if max(Fid)==5 || max(Fid)==6
    basal_nodes=find(labelp==4);
elseif max(Fid)==18
    %basal_nodes=find(labelp==9 | labelp==10 | labelp==13 | labelp==14 );
   basal_nodes= find(a2b >=Threshold_val);
end
% Laplace
boundaryIds=[C_apex_points_id; basal_nodes];
boundaryVal=[zeros(size(C_apex_points_id)); ones(size(basal_nodes))];

a2b_cobi = solveLaplace(L, boundaryIds, boundaryVal, 1e-8, 5000);
abLaplaceGrad_c = normalizedGradField(GT, a2b_cobi , 1e-8, false, v, elem);


%%
    e0_elem=cross(tmGrad_c,abLaplaceGrad_c);
     warning('off','all')

    %%Septum
    indx_sR_post_relative=find(ismember(Sept_points,indx_post_row));
    indx_cells_septum=ismember(elem,elem(idx_sept_elems,:),'rows');
    indx_cells_septum2=and(index_heart_elem,indx_cells_septum);
    dSeptPost = solveTrajectDist(GSept, e0_elem(indx_cells_septum2,:),  indx_sR_post_relative, zeros(size( indx_sR_post_relative)), 1e-8, 3000);
    
    indx_sR_ant_relative=find(ismember(Sept_points,indx_ant_points));
    dSeptAnt = solveTrajectDist(GSept, -e0_elem(indx_cells_septum2,:),  indx_sR_ant_relative, zeros(size( indx_sR_ant_relative)), 1e-9, 3000);
    % write_vtk_rbm('Testdd.vtk',length(Sept_points),v(Sept_points,:),length(Sept_elem),wSeptF,  dSeptAnt); %write mesh

    rtTrajectDistSept = dSeptPost./(dSeptAnt+dSeptPost);
     rtTrajectDistSept( rtTrajectDistSept>1)=1;
     rtTrajectDistSept( rtTrajectDistSept<0)=0;
    %fix points of boundaries 

    rtSept = 2/3 + 1/3 * (1-rtTrajectDistSept);

    %%Free
    indx_fR_post_relative=find(ismember(Free_points,indx_post_row_free));
    indx_cells_free=ismember(elem,elem(~idx_sept_elems,:),'rows');
    indx_cells_free2=and(index_heart_elem,indx_cells_free);%filter base elements
    dFreePost = solveTrajectDist(GFree, e0_elem(indx_cells_free2,:),  indx_fR_post_relative, zeros(size( indx_fR_post_relative)), 1e-8, 3000);

    indx_fR_ant_relative=find(ismember(Free_points,indx_ant_points_free));
    dFreeAnt = solveTrajectDist(GFree, -e0_elem(indx_cells_free2,:),  indx_fR_ant_relative, zeros(size( indx_fR_ant_relative)), 1e-8, 3000);

    rtTrajectDistFree = dFreePost./(dFreeAnt+dFreePost);
    rtTrajectDistFree( rtTrajectDistFree>1)=1;
    rtTrajectDistFree( rtTrajectDistFree<0)=0;
    rtFree = 2/3 * rtTrajectDistFree;
    r = NaN(size(v,1),1);
    r(Free_points) = rtFree;
    r(Sept_points) = rtSept;   
    rsin=sin(2*pi*r);
    rcos=cos(2*pi*r);   
    r3=r;
    r3(isnan(r))=-1; %random value





 %%    a2b coordinate
 %we use an adapted function from the original cobiveco code
 %Input preparation
 index_heart_point=unique(elem(index_heart_elem,:));
 oo.m2.vol.points=v(index_heart_point,:);
 oo.m2.vol.cells=int32(sortconnectivities(elem(index_heart_elem,:)));
 oo.m2.vol.cellTypes=uint8(ones(length(elem(index_heart_elem,:)),1)*10);
 oo.m2.rtSin=rsin(index_heart_point);
 oo.m2.rtCos=rcos(index_heart_point);
 oo.m2.r=r3(index_heart_point);
 P1 = double(oo.m2.vol.points);
 C1 = double(oo.m2.vol.cells);
 oo.m2.sur = vtkDataSetSurfaceFilter((oo.m2.vol));
 oo.m2.surToVol = vtkMapPointIds(oo.m2.vol, oo.m2.sur);

 %laplace values are modified in case we use a cut geometry
if max(Fid)==5 || max(Fid)==6
     L_ab=L;
     GT_ab=GT;
     oo.m2.abLaplace=a2b_cobi(index_heart_point);
     %basal nodes
     indx_row=(label==4);
     indx_nodes_surf_base=unique(oo.m2.sur.cells(indx_row,:));

elseif max(Fid)>=18 
     L_ab=cotmatrix(oo.m2.vol.points, oo.m2.vol.cells);
     GT_ab= grad(double(oo.m2.vol.points),double( oo.m2.vol.cells));
    %new basal nodes
    [~,indx_nodes_surf_base]=setdiff(oo.m2.sur.points(:,:),v(unique(face),:),"rows");
     oo.m2.abLaplace=a2b_cobi(index_heart_point);

end
 oo.m2.L = L_ab;
 oo.m2.G = GT_ab;
 oo.m2.massMat = massmatrix(P1, C1, 'voronoi');
 oo.m2.tm=tm_cobi(index_heart_point);
 oo.m2.meanEdgLen=mean(vtkEdgeLengths(oo.m2.vol));
 class=zeros(length(oo.m2.sur.points),1);  
 class(indx_nodes_surf_base)=1;%basal label
 oo.m2.sur.pointData.class= class;
 oo.m1=oo.m2;

 a2b_cobi_coord=computeApicobasal_adapted(oo);

warning('on','all')
 if max(Fid==18)
    a2b_cobi_coord2=a2b_cobi_coord;
    a2b_cobi_coord=-1*ones(size(r));
    a2b_cobi_coord(index_heart_point)=a2b_cobi_coord2;
 end
  
  a2b_cobi_coord_final=a2b_cobi_coord;
  a2b_cobi_coord_final(isnan(r))=NaN; %random value
  %%   
  %f=[ones(length(elem),1)*4   (elem)-1];
  %write_vtk_rbm('C_test.vtk',length(v),v,length(f),f,tm_cobi, a2b_cobi_coord,r3)%write mesh


 %function
    function  Sorted_Con=sortconnectivities(connectivities_orig)
            %reorder mesh connectivities
                [a1, ~, a2] = unique(connectivities_orig);
                Anew=1:length(a1);
                Surfaces0 = Anew(a2);
                Sorted_Con= reshape(Surfaces0, size(connectivities_orig));
    end

end