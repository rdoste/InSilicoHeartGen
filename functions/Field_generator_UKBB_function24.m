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


function Field_generator_UKBB_function24(Fiber_info,meshformat,pericardiumlevel, epiendo, epiendoRV, case_number)
%add original projected cobiveco coodinates
directory=pwd;

%Paths=====================================================================================
directoryResults=pwd;


%% Input variables=================================================================================================
  %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     name='coarse.vtu'; %name of the tetrahedral mesh (in case largenumber==1, this is the coarse mesh)
     labelmesh='labels_final.vtk'; %mesh with labels
     name_large='fine.vtu'; %name of the hexahedral or the final tetrahedral mesh (only if largenumber is activated)
     largenumber=1;  %0--> no need interpolation for a large number of points
                     %1--> interpolation requiered after Elmer                                                
      %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    %% First Simulation   (only for calculate transmural gradient and determine RV and LV)
    %open vtk=================================================================================================

     [~,pto,~,car,Fid]=read_labels2(labelmesh);           
     [MeshCoarse]=vtkRead(name);
     node=double(MeshCoarse.points);
     elem=double(MeshCoarse.cells);
     face=volface(elem);

    %labelling volumetric mesh (define boundaries for Heat Equation)
    disp('Labelling and mesh preparation');

             centroid=double(meshcentroid(node,face));
             centroid_lab=meshcentroid(pto,car+1);
             TR=delaunayTriangulation(centroid_lab);
             TR_mesh=triangulation(double(elem),double(node));

             [N,~]=nearestNeighbor(TR,centroid);
             label=Fid(N);
             %Force at least one face as apex
              %cut geometry
           if meshformat=="cut"
              [~,Apex_ID]= min(sqrt(sum((centroid - reshape(centroid_lab(Fid==5,:)',1,3,[])).^2,2))); % Distance from apex centroid
              label(Apex_ID)=5;
           else
              %biventricular
              [~,ApexLV_ID]= min(sqrt(sum((centroid - reshape(centroid_lab(Fid==12,:)',1,3,[])).^2,2))); % Distance from apex centroid
              label(ApexLV_ID)=12;
              [~,ApexRV_ID]= min(sqrt(sum((centroid - reshape(centroid_lab(Fid==18,:)',1,3,[])).^2,2))); % Distance from apex centroid
              label(ApexRV_ID)=18;
           end



             centroid_tetra=meshcentroid(node,elem);
             [N_tetra,~]=nearestNeighbor(TR,centroid_tetra);
             label_tetra=Fid(N_tetra);
             label_tetra2=ones(size(label_tetra)); % we force only one physical body
             nnode=length(node);
             nelem=length(elem);
             carelem=size(elem,2); 
        
          %label projection to nodes
            [N2,dist]=nearestNeighbor(TR,node);         
            labelp=Fid(N2);

            %%make sure to add apex points

              %%make sure to add apex points
         if meshformat=="cut"
            labelp(face(label==5,:))=5;
            labelp(face(label==4,:))=4;
         else
            labelp(face(label==12,:))=12;
            labelp(face(label==18,:))=18;
         end
            write_vtk_surf('labelspoints.vtk',node,face,labelp);


        %%
        disp('Solving Laplace equations...')
            L=cotmatrix(node, elem);
            GT= grad(double(node), double(elem));

%%

        %Transmural****************************************************************************************************************************

           disp('Computing Transmural gradients...')
            id1=intersect(find(ismember(labelp,[1 5 12 19 20 23 24])),unique(face));
            id2=intersect(find(ismember(labelp,[3 6 9 14 18])),unique(face));
            id3=intersect(find(ismember(labelp,[2 10 13])),unique(face));
            ids = [id1;id2;id3];
            val = [ zeros(size(id1)); ones(size(id2)); -2.*ones(size(id3))];
            FieldsC.Tphi = solveLaplace(L, ids, val, 1e-9, 5000);
            transmural_grad = normalizedGradField(GT, FieldsC.Tphi,1e-10, false, node, elem);
            FieldsC.DTphi=elemGrad2pointGrad(transmural_grad,TR_mesh); % projection to mesh nodes
            
            % elem_w=[ones(length(elem),1)*4   (elem)-1];
            % write_vtk_rbm('E_trans.vtk',length(node),node,length(elem),elem_w,FieldsC.Tphi,FieldsC.DTphi); %write mesh
        %Transmural_bi ************************************************************************************************************************
            id1=intersect(find(ismember(labelp,[1 5 12 19 20 23 24])),unique(face));
            id2=intersect(find(ismember(labelp,[3 6 9 14 18])),unique(face));
            id3=intersect(find(ismember(labelp,[2 10 13])),unique(face));
            ids = [id1;id2;id3];
            val = [ zeros(size(id1)); ones(size(id2)); ones(size(id3))];
            FieldsC.Tphi_bi = solveLaplace(L, ids, val, 1e-9, 5000);
            transmural_bi_grad = normalizedGradField(GT, FieldsC.Tphi_bi,1e-10, false, node, elem);
            FieldsC.DTphi_bi=elemGrad2pointGrad(transmural_bi_grad,TR_mesh); % projection to mesh nodes
             % write_vtk_rbm('E_transbi.vtk',length(node),node,length(elem),elem_w,FieldsC.Tphi_bi,FieldsC.DTphi_bi); %write mesh
         %Transventricular_Cobi ************************************************************************************************************************

            id1=intersect(find(ismember(labelp,[3 6 9 19 14 24 18])),unique(face));
            id2=intersect(find(ismember(labelp,[2 10 20 13 23])),unique(face));
            ids = [id1;id2];
            val = [ ones(size(id1)); zeros(size(id2))];
            FieldsC.Tv_cobi = solveLaplace(L, ids, val, 1e-9, 5000);
            transventricular_cobi_grad = normalizedGradField(GT, FieldsC.Tv_cobi,1e-10, false, node, elem);
            FieldsC.DTv_cobi=elemGrad2pointGrad(transventricular_cobi_grad,TR_mesh); % projection to mesh nodes
            % write_vtk_rbm('E_transventr.vtk',length(node),node,length(elem),elem_w,FieldsC.Tv_cobi,FieldsC.DTv_cobi); %write mesh


         %Transmural RV****************************************************************************************************************************

            id1=intersect(find(ismember(labelp,[1 5 6 12 19 20 23 24])),unique(face));
            id2=intersect(find(ismember(labelp,[3 9 14 18])),unique(face));
            id3=intersect(find(ismember(labelp,[2 10 13])),unique(face));
            ids = [id1;id2;id3];
            val = [ zeros(size(id1)); ones(size(id2)); ones(size(id3))];
            FieldsC.Tphi3 = solveLaplace(L, ids, val, 1e-9, 5000);
            transmural3_grad = normalizedGradField(GT, FieldsC.Tphi3,1e-10, false, node, elem);
            FieldsC.DTphi3=elemGrad2pointGrad(transmural3_grad,TR_mesh); % projection to mesh nodes
            % write_vtk_rbm('E_trans3.vtk',length(node),node,length(elem),elem_w,FieldsC.Tphi3,FieldsC.DTphi3); %write mesh

           % Read transmural information and Ventricle definition                           


            N_points=length(node);
            v=node;
            f=[ones(length(elem),1)*4   (elem)-1];
            N_faces=length(f);
            f2=f; %save original faces
            v2=v;    %save original (future aha calculation)          


            %% create Fid again with old valve values %%%%%%%%%%%%%%%


            %%make sure to add apex points
            FieldsC.DTphi=double(normr(FieldsC.DTphi));

        %%%%%% We define LV or RV tags (epi-endo)%%%%%%%%%%%%%%%%%%%%%%%%%
          Ventricle=single(FieldsC.Tphi); 
          Ventricle(Ventricle<0 & Ventricle>-2)=-1; % mio LV =-1
          Ventricle(Ventricle>0 & Ventricle<1)=2;  % mio RV =2 
          gradientdot=dot(FieldsC.DTphi,FieldsC.DTphi_bi,2);%the different direction of the gradients will determine LV vs RV
          Ventricle(FieldsC.Tphi==0 & gradientdot>0)=3;
          Ventricle(FieldsC.Tphi==0 & gradientdot<0)=-3;

          FieldsC.Ventricle=Ventricle;
%% new msh creation
              body=FieldsC.Tphi(elem)./abs(FieldsC.Tphi(elem));
              body1=mode(body,2);


       %RV mesh
        elemRV=elem(body1>0,:);
        indp_RV=unique(elemRV);
        L_RV=cotmatrix(node(indp_RV,:), sortconnectivities(elemRV));
        GT_RV= grad(double(node(indp_RV,:)), double(sortconnectivities(elemRV)));
        TR_mesh_RV=triangulation(sortconnectivities(elemRV),node(indp_RV,:));
      %LV mesh
        elemLV=elem(body1<0,:);
        indp_LV=unique(elemLV);
        L_LV=cotmatrix(node(indp_LV,:), sortconnectivities(elemLV));
        GT_LV= grad(double(node(indp_LV,:)), double(sortconnectivities(elemLV)));
        TR_mesh_LV=triangulation(sortconnectivities(elemLV),node(indp_LV,:));

      %check that the apex is in the LV
      if meshformat~="cut"
           id1=intersect(find(ismember(labelp,12)),unique(face));
           ind1_rel=find(ismember(indp_LV,id1),1);
           if isempty(ind1_rel) %if apex is in the RV, recalculate it
               label(label==12)=1;
               centerLV=mean(centroid(label==2,:));
               centerLVring=mean(centroid(label==10 | label==13,:));
               LVaxis=centerLV-centerLVring;
               heightLV = double(node*LVaxis');
               [~,ApexLV_id_LV]=max(heightLV);
               labelp(labelp==12)=1;
               labelp(ApexLV_id_LV)=12;
               [~,ApexLV_id_LV_cells]=ismember(ApexLV_id_LV,face);
               label(ApexLV_id_LV_cells)=12;
           end
      end
        

      

        %% only for closed meshes
        Fid2=Fid;
        Fid(Fid>=19)=Fid(Fid>=19)-10;
        label_closed=label;
        label(label>=19)= label(label>=19)-10;
        labelp(labelp>=19)= labelp(labelp>=19)-10;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

%%
     %%Longitudinal_bi***********************************************************************************************************************
        disp('Computing longitudinal gradients...')
        id1=intersect(find(ismember(labelp,[5 12])),unique(face));
        id2=intersect(find(ismember(labelp,[4 9 14 10 13])),unique(face));
        ids = [id1;id2];
        val = [ ones(size(id1)); zeros(size(id2))];
        FieldsC.Tpsi_bi = solveLaplace(L, ids, val, 1e-9, 5000);
        longitudinal_bi_grad = normalizedGradField(GT, FieldsC.Tpsi_bi,1e-10, false, node, elem);
        FieldsC.DTpsi_bi=elemGrad2pointGrad(longitudinal_bi_grad,TR_mesh); % projection to mesh nodes
        % write_vtk_rbm('E_long.vtk',length(node),node,length(elem),elem_w,FieldsC.Tpsi_bi,FieldsC.DTpsi_bi); %write mesh

   if meshformat~="cut"
      %Long_LR*******************************************************************************************************************************
        id1=intersect(find(ismember(labelp,18)),unique(face));
        [ind1_rel]=find(ismember(indp_RV,id1));
        id2=intersect(find(ismember(labelp,[9 14])),unique(face));
        [ind2_rel]=find(ismember(indp_RV,id2));
        ids = [ind1_rel;ind2_rel];
        val = [ ones(size(ind1_rel)); zeros(size(ind2_rel))];
        FieldsC.TpsiR0 = solveLaplace(L_RV, ids, val, 1e-9, 5000);
        longitudinal_bi_grad = normalizedGradField(GT_RV, FieldsC.TpsiR0,1e-10, false, node,elemRV);
        FieldsC.DTpsiR0=elemGrad2pointGrad(longitudinal_bi_grad,TR_mesh_RV); % projection to mesh nodes
        % elem_w_RV=[ones(length(elem(body1>0,:)),1)*4   (sortconnectivities(elemRV))-1];
        % write_vtk_rbm('E_longRV.vtk',length(node(indp_RV,:)),node(indp_RV,:),length(elemRV),elem_w_RV,FieldsC.TpsiR0,FieldsC.DTpsiR0); %write mesh



      %Long_LV*******************************************************************************************************************************
        id1=intersect(find(ismember(labelp,12)),unique(face));
        [ind1_rel]=find(ismember(indp_LV,id1));
        id2=intersect(find(ismember(labelp,[10 13])),unique(face));
        [ind2_rel]=find(ismember(indp_LV,id2));
        ids = [ind1_rel;ind2_rel];
        val = [ ones(size(ind1_rel)); zeros(size(ind2_rel))];
        FieldsC.TpsiL0 = solveLaplace(L_LV, ids, val, 1e-9, 5000);
        longitudinal_bi_grad = normalizedGradField(GT_LV, FieldsC.TpsiL0,1e-10, false, node,elemLV);
        FieldsC.DTpsiL0=elemGrad2pointGrad(longitudinal_bi_grad,TR_mesh_LV); % projection to mesh nodes
        % elem_w_LV=[ones(length(elemLV),1)*4   (sortconnectivities(elemLV))-1];
        % write_vtk_rbm('E_longLV.vtk',length(node(indp_LV,:)),node(indp_LV,:),length(elemLV),elem_w_LV,FieldsC.TpsiL0,FieldsC.DTpsiL0); %write mesh

      % 2valveRV*******************************************************************************************************************************
        id1=intersect(find(ismember(labelp,18)),unique(face));
        [ind1_rel]=find(ismember(indp_RV,id1));
        id2=intersect(find(ismember(labelp,14)),unique(face));
        [ind2_rel]=find(ismember(indp_RV,id2));
        id3=intersect(find(ismember(labelp,9)),unique(face));
        [ind3_rel]=find(ismember(indp_RV,id3));
        ids = [ind1_rel;ind2_rel;ind3_rel];
        val = [ ones(size(ind1_rel)); ones(size(ind2_rel)); zeros(size(ind3_rel))];
        FieldsC.Tval0= solveLaplace(L_RV, ids, val, 1e-9, 5000);  
        % elem_w_RV=[ones(length(elem(body1>0,:)),1)*4   (sortconnectivities(elemRV))-1];
        % write_vtk_rbm('E_2vRV.vtk',length(node(indp_RV,:)),node(indp_RV,:),length(elemRV),elem_w_RV,FieldsC.Tval0); %write mesh
     %2valveLV*******************************************************************************************************************************

        id1=intersect(find(ismember(labelp,12)),unique(face));
        [ind1_rel]=find(ismember(indp_LV,id1));
        id2=intersect(find(ismember(labelp,13)),unique(face));
        [ind2_rel]=find(ismember(indp_LV,id2));
        id3=intersect(find(ismember(labelp,10)),unique(face));
        [ind3_rel]=find(ismember(indp_LV,id3));
        ids = [ind1_rel;ind2_rel;ind3_rel];
        val = [ ones(size(ind1_rel)); ones(size(ind2_rel)); zeros(size(ind3_rel))];
        FieldsC.Tval_l0= solveLaplace(L_LV, ids, val, 1e-9, 5000);  
        % elem_w_LV=[ones(length(elemLV),1)*4   (sortconnectivities(elemLV))-1];
        % write_vtk_rbm('E_2vLV.vtk',length(node(indp_LV,:)),node(indp_LV,:),length(elemLV),elem_w_LV,FieldsC.Tval_l0); %write mesh



   end

    %% location of the septum surface creating two bodies (one for each
   %ventricle)

   %Septum****************************************************************************************************************************
     disp('Computing septal gradients...')
           
              [face2,labelf2,label_tetra_msh]=msh_septum_creation(v2,elem,face,labelp,FieldsC.Tphi,label,body1,label_tetra,0,20,22);
                %write_vtk_surf('labelsS.vtk',v2,face2,labelf2);

       %labels to points
        centroid_lab2=meshcentroid(v2,face2);
        warning('off','all')
        TR2=delaunayTriangulation(centroid_lab2);
        warning('on','all')
        [N2,dist]=nearestNeighbor(TR2,v2);         
        labelp2=labelf2(N2);
        labelp2(unique(face2(labelf2==20,:)))=20;%force that all point of the septal surface  have the same value (including epicardial)

        id1=intersect(find(ismember(labelp2,[3 6 2 18])),unique(face2));
        id2=intersect(find(ismember(labelp2,20)),unique(face2));
        ids = [id1;id2];
        val = [ zeros(size(id1)); ones(size(id2))];
        FieldsC.T_sept = solveLaplace(L, ids, val, 1e-9, 5000);
        septum_grad = normalizedGradField(GT,  FieldsC.T_sept,1e-10, false, node, elem);
        FieldsC.DT_sept=elemGrad2pointGrad( septum_grad,TR_mesh); % projection to mesh nodes
        % write_vtk_rbm('E_septum.vtk',length(node),node,length(elem),elem_w, FieldsC.T_sept, FieldsC.DT_sept); %write mesh
    
   

   %Septum Cobiveco***************************************************************************************************************************

              body01=FieldsC.Tv_cobi(elem);
              body11=mode(body01,2);
              [face3,labelf3,label_tetra_msh3]=msh_septum_creation2(v2,elem,face,labelp2,FieldsC.Tphi,label_closed,body11,label_tetra,0.5,30,32);
             % write_vtk_surf('labelsS2.vtk',v2,face3,labelf3);

    %Septum****************************************************************************************************************************

            %labels to points
            centroid_lab3=meshcentroid(v2,face3);
            warning('off','all')
            TR3=delaunayTriangulation(centroid_lab3);
            warning('on','all')
            [N2,dist]=nearestNeighbor(TR3,v2);         
            labelp3=labelf3(N2);
            labelp3(unique(face3(labelf3==30,:)))=30;%force that all point of the septal surface  have the same value (including epicardial)

            id1=intersect(find(ismember(labelp3,[1 5 12])),unique(face3));
            id2=intersect(find(ismember(labelp3,30)),unique(face3));
            ids = [id1;id2];
            val = [ zeros(size(id1)); ones(size(id2))];
            FieldsC.T_sept_cobi = solveLaplace(L, ids, val, 1e-9, 5000);
            septum_cobi_grad = normalizedGradField(GT,  FieldsC.T_sept_cobi,1e-10, false, node, elem);
            FieldsC.DT_sept_cobi=elemGrad2pointGrad( septum_cobi_grad,TR_mesh); % projection to mesh nodes
            % write_vtk_rbm('E_septum_cobi.vtk',length(node),node,length(elem),elem_w, FieldsC.T_sept_cobi, FieldsC.DT_sept_cobi); %write mesh


            id1=intersect(find(ismember(labelp3,[3 6 2 18])),unique(face3));
            id2=intersect(find(ismember(labelp3,30)),unique(face3));
            ids = [id1;id2];
            val = [ zeros(size(id1)); ones(size(id2))];
            FieldsC.T_sept_uvc = solveLaplace(L, ids, val, 1e-9, 5000);
            septum_uvc_grad= normalizedGradField(GT,  FieldsC.T_sept_uvc,1e-10, false, node, elem);
            FieldsC.DT_sept_uvc=elemGrad2pointGrad( septum_uvc_grad,TR_mesh); % projection to mesh nodes
            % write_vtk_rbm('E_septum_uvc.vtk',length(node),node,length(elem),elem_w, FieldsC.T_sept_uvc, FieldsC.DT_sept_uvc); %write mesh


      %Transmural_cobi****************************************************************************************************************************

            id1=intersect(find(ismember(labelp3,[1 5 12 19 20 23 24 30])),unique(face3));
            id2=intersect(find(ismember(labelp3,[2 3 6 9 10 13 14 18])),unique(face3));
            ids = [id1;id2];
            val = [ zeros(size(id1)); ones(size(id2))];
            FieldsC.Tphi_cobi = solveLaplace(L, ids, val, 1e-9, 5000);
            Transmural_cobi_grad = normalizedGradField(GT,  FieldsC.Tphi_cobi,1e-10, false, node, elem);
            FieldsC.DTphi_cobi=elemGrad2pointGrad( Transmural_cobi_grad,TR_mesh); % projection to mesh nodes
            % write_vtk_rbm('E_Transmural_cobi.vtk',length(node),node,length(elem),elem_w, FieldsC.Tphi_cobi, FieldsC.DTphi_cobi); %write mesh


         %% Fiber Generation 
            
%%
      
%Transfer some of the data calculated in 1 ventricle to the whole heart
   indR=find(FieldsC.Ventricle==1 | FieldsC.Ventricle==2 | FieldsC.Ventricle==3); %RV
   indL=find(FieldsC.Ventricle==-1 | FieldsC.Ventricle==-2 | FieldsC.Ventricle==-3); %LV

   [indR2,indR2_pos]=intersect(indp_RV,indR);
   [indL2,indL2_pos]=intersect(indp_LV,indL);

if meshformat~="cut"
   FieldsC.Tval=zeros(size(FieldsC.Tphi));      FieldsC.Tval(indR2)=FieldsC.Tval0(indR2_pos);            FieldsC.Tval(indL2)=NaN;
   FieldsC.Tval_l=zeros(size(FieldsC.Tphi));    FieldsC.Tval_l(indL2)=FieldsC.Tval_l0(indL2_pos);        FieldsC.Tval_l(indR2)=NaN;     
   FieldsC.TpsiR=zeros(size(FieldsC.Tphi));     FieldsC.TpsiR(indR2)=FieldsC.TpsiR0(indR2_pos);          FieldsC.TpsiR(indL2)=NaN;
   FieldsC.TpsiL=zeros(size(FieldsC.Tphi));     FieldsC.TpsiL(indL2)=FieldsC.TpsiL0(indL2_pos);          FieldsC.TpsiL(indR2)=NaN;     
   FieldsC.DTpsiR=zeros(size(FieldsC.DTphi));   FieldsC.DTpsiR(indR2,:)=FieldsC.DTpsiR0(indR2_pos,:);    FieldsC.DTpsiR(indL2,:)=repmat([0 0 0],length(indL2),1);
   FieldsC.DTpsiL=zeros(size(FieldsC.DTphi));   FieldsC.DTpsiL(indL2,:)=FieldsC.DTpsiL0(indL2_pos,:);    FieldsC.DTpsiL(indR2,:)=repmat([0 0 0],length(indR2),1);
end

    %normalize gradients and change sign in fluxes
    
if meshformat~="cut"
    % FieldsC.DTpsiL=double(normr(FieldsC.DTpsiL));
    % FieldsC.DTpsiR=double(normr(FieldsC.DTpsiR));
    FieldsC.Tval(FieldsC.Tval<0)=0;
    FieldsC.Tval_l(FieldsC.Tval_l<0)=0;
end
    FieldsC.DTphi=double(normr(FieldsC.DTphi));
    FieldsC.Tphi3=min(max(FieldsC.Tphi3,0),1); %filter some problematic values out of range in the boundary

    Fields=FieldsC;

   %% Interpolation for large numbers
                 if largenumber==1
                     disp('Interpolation to finer mesh');

                       [MeshFine]=vtkRead(name_large);
                      
                        node_vox=double(MeshFine.points);                 
                        elem_vox=MeshFine.cells;
                      
                        % Laplace solutions interpolation

                        T=triangulation(double(f2(:,2:5)+1),double(v2));
                       
                        [tri,bar]=pointLocation(T,node_vox(:,:));
                        
                        %substitute points outside the coarse mesh (NaN) by
                        %nearestNeighbor points values
                        [row, ~] = find(isnan(tri));
                        pointN=(nearestNeighbor(T,node_vox(row,:)));                      
                        [~,n2]=ismember(pointN,f2(:,2:end)+1);
                        [row3, ~] = ind2sub(size(f2(:,2:end)), n2);
                        tri(row)=row3;                                               
                        bar(row,:)=cartesianToBarycentric(T,row3,v2(pointN,:));
                    
                        points_Tetra=f2(tri,2:end)+1;

                                              
                        % Interpolation
                       
                        TphipointsTetra= FieldsC.Tphi(points_Tetra);
                        FieldsF.Tphi=dot(bar,TphipointsTetra,2);
                        FieldsF.Tphi(row)=round( FieldsF.Tphi(row)); % remove interpolation artifacts in border if NaN

                        TphibipointsTetra=FieldsC.Tphi_bi(points_Tetra);
                        FieldsF.Tphi_bi=dot(bar,TphibipointsTetra,2);
                        FieldsF.Tphi_bi(row)=round(FieldsF.Tphi_bi(row)); % remove interpolation artifacts in border if NaN

                        TphicobipointsTetra=FieldsC.Tphi_cobi(points_Tetra);
                        FieldsF.Tphi_cobi=dot(bar,TphicobipointsTetra,2);
                        FieldsF.Tphi_cobi(row)=round(FieldsF.Tphi_cobi(row)); % remove interpolation artifacts in border if NaN
                
                        TpsibipointsTetra=FieldsC.Tpsi_bi(points_Tetra);
                        FieldsF.Tpsi_bi=dot(bar,TpsibipointsTetra,2);

                        TseptpointsTetra=FieldsC.T_sept(points_Tetra);
                        FieldsF.T_sept=dot(bar,TseptpointsTetra,2);

                        Tsept_cobipointsTetra=FieldsC.T_sept_cobi(points_Tetra);
                        FieldsF.T_sept_cobi=dot(bar,Tsept_cobipointsTetra,2);

                        Tsept_uvcpointsTetra=FieldsC.T_sept_uvc(points_Tetra);
                        FieldsF.T_sept_uvc=dot(bar,Tsept_uvcpointsTetra,2);

                        Tv_cobipointsTetra=FieldsC.Tv_cobi(points_Tetra);
                        FieldsF.Tv_cobi=dot(bar,Tv_cobipointsTetra,2);
                
                  if meshformat~="cut"
                  
                        TpsiRVpointsTetra=FieldsC.TpsiR(points_Tetra);
                        % FieldsF.TpsiR=dot(bar,TpsiRVpointsTetra,2);
                        FieldsF.TpsiR=baryInterpIgnoreNaN(bar, TpsiRVpointsTetra); %interpolation between ventricles
                        FieldsF.TpsiR(isnan(FieldsF.TpsiR))=-1;
                        
                        TpsiLVpointsTetra=FieldsC.TpsiL(points_Tetra);
                        % FieldsF.TpsiL=dot(bar,TpsiLVpointsTetra,2);
                        FieldsF.TpsiL=baryInterpIgnoreNaN(bar,TpsiLVpointsTetra); %interpolation between ventricles
                        FieldsF.TpsiL(isnan(FieldsF.TpsiL))=-1;


                        TvalpointsTetra=FieldsC.Tval(points_Tetra);
                        % FieldsF.Tval=dot(bar,TvalpointsTetra,2);
                        FieldsF.Tval=baryInterpIgnoreNaN(bar,TvalpointsTetra); %interpolation between ventricles
                        FieldsF.Tval(isnan(FieldsF.Tval))=-1;

                        
                        TvallpointsTetra=FieldsC.Tval_l(points_Tetra);
                        % FieldsF.Tval_l=dot(bar,TvallpointsTetra,2);
                        FieldsF.Tval_l=baryInterpIgnoreNaN(bar, TvallpointsTetra); %interpolation between ventricles
                        FieldsF.Tval_l(isnan(FieldsF.Tval_l))=-1;


                  end

                        Tphi3pointsTetra=FieldsC.Tphi3(points_Tetra);
                        FieldsF.Tphi3=dot(bar,Tphi3pointsTetra,2);
                        FieldsF.Tphi3(row)=round(FieldsF.Tphi3(row)); % remove interpolation artifacts in border if NaN

                        %  Gradient interpolation
                       
                         FieldsF.DTphi=gradientinterpol(FieldsC.DTphi,length(node_vox),points_Tetra,bar);
                         FieldsF.DTphi_bi=gradientinterpol(FieldsC.DTphi_bi,length(node_vox),points_Tetra,bar);
                         FieldsF.DTpsi_bi=gradientinterpol(FieldsC.DTpsi_bi,length(node_vox),points_Tetra,bar);
                  if meshformat~="cut"       
                         FieldsF.DTpsiR=gradientinterpol(FieldsC.DTpsiR,length(node_vox),points_Tetra,bar);
                         FieldsF.DTpsiL=gradientinterpol(FieldsC.DTpsiL,length(node_vox),points_Tetra,bar);
                  end
                        %change points of calculations
                        v=node_vox;
                        N_points=length(node_vox);
                        f=[ones(length(elem_vox),1)*4,elem_vox-1];  %remove in the future
                        N_faces=length(elem_vox);
                        name=name_large;
                
                %Ventricle Tag interpolation
                        
                          % Ventricle_aux=int8(Ventricle(points_Tetra));
                          % Ventricle_raw=round(dot(bar,double(Ventricle_aux),2));%know lV vs RV 
                          % %rounded value
                          % Epiendo3=Ventricle_aux;%auxiliar
                          % %positive values
                          % Epiendo3(Ventricle_aux==-2)=1;
                          % Epiendo3(Ventricle_aux==-1)=2;
                          % Epiendo3(Ventricle_aux==-3)=3;
                          % %Ventricle values calculation
                          % Ventricle=round(dot(bar,double(Epiendo3),2));
                          % Ventricle(Ventricle_raw<0)=Ventricle(Ventricle_raw<0).*-1;            
                          % Epiendo3=int8(Ventricle);%auxiliar
                          % Ventricle(Epiendo3==-1)=-2;
                          % Ventricle(Epiendo3==-2)=-1;
                          % FieldsF.Ventricle=Ventricle;


                          face_fine=volface(f(:,2:5)+1);
                          surf_nodes_ind=sort(unique(face_fine));%index of surface nodes
                          surf_nodes = zeros(length(v),1); surf_nodes(surf_nodes_ind) = 1; %in boolean
                          FieldsF.Ventricle=single(FieldsF.Tphi); 
                          FieldsF.Ventricle(FieldsF.Ventricle<0 & FieldsF.Ventricle>-2)=-1; % myo LV =-1
                          FieldsF.Ventricle(FieldsF.Ventricle>=0 & FieldsF.Ventricle<1)=2;  % myo RV =2 
                          gradientdot=dot(FieldsF.DTphi,FieldsF.DTphi_bi,2);%the different direction of the gradients will determine LV vs RV
                          FieldsF.Ventricle(surf_nodes & FieldsF.Tphi_bi<0.1 & gradientdot>0)=3;
                          FieldsF.Ventricle(surf_nodes & FieldsF.Tphi_bi<0.1  & gradientdot<0)=-3;
                          FieldsF.Ventricle(surf_nodes & FieldsF.Tphi>0.9  & gradientdot>0)=1;
                          FieldsF.Ventricle(surf_nodes & FieldsF.Tphi<-1.9  & gradientdot<0)=-2;
                   clear -regexp pointsTetra$; 

                   Fields=FieldsF;
                        
                 end
                
  %% Axis
  disp('Fiber Calculation');
    % e1 longitudinal+vertical(apex 2 mitral and triculpid valves)
    if meshformat=="cut"
        e1=normr(Fields.DTpsi_bi);
        e1=normr(e1);
    else
        Fields.DTpsi_result=normr(Fields.DTpsiL);  %2ValvesDirection+PVDirection
        indR=find(Fields.Ventricle==1 | Fields.Ventricle==2 | Fields.Ventricle==3); %RV
        Fields.DTpsi_result(indR,:)=normr(Fields.DTpsiR(indR,:));
        e1=normr(Fields.DTpsi_result);
        e1=normr(e1);

    end

    %%%%%% FieldsC.DTphi_final is an interpolation between the two transmural directions
    %try bislerp (to do)
        Fields.DTphi(Fields.Ventricle<0,:)=Fields.DTphi(Fields.Ventricle<0,:).*-1;
        Fields.DTphi_final=normr(Fields.Tphi_bi.*Fields.DTphi+(1-Fields.Tphi_bi).*Fields.DTphi_bi);
        Fields.DTphi_final(Fields.Ventricle<0,:)=Fields.DTphi_final(Fields.Ventricle<0,:).*-1;
        proj=Fields.DTphi_final-bsxfun(@times,dot(e1',Fields.DTphi_final')',e1);
    %%%%%%%%%%%%

        e2=normr(proj);
        e2=normr(e2);
        e0=cross(e1',e2')';
        e0=normr(e0);

    %Define transmurality (0-epi   1-endo)
    d=abs(abs(Fields.Tphi)-1); %invert. d=0 endo d=1 epi
    d2=abs(abs(Fields.Tphi)/2-1); %for LV
    d_bi=Fields.Tphi_bi;
if meshformat=="cut"
    Fields.Tval_l=ones(size(Fields.Tphi_bi));  %force Field.Tval to be ones, not affecting the angle computation
    Fields.Tval=ones(size(Fields.Tphi_bi));
end
%%  alpha angle 

    alpha_wall=zeros(size(Fields.Ventricle));
    alpha_wall_epi=zeros(size(Fields.Ventricle));
    alpha_wall_endo=zeros(size(Fields.Ventricle));     
    
    %Epicardium orientation             
          alpha_wall(Fields.Ventricle==-3)=deg2rad(Fiber_info.AEPILV).*Fields.Tval_l(Fields.Ventricle == -3) +deg2rad(Fiber_info.AOTEPILV) .*(1-Fields.Tval_l(Fields.Ventricle==-3));   %OT interpolation
          alpha_wall(Fields.Ventricle==3 )=deg2rad(Fiber_info.AEPIRV).*Fields.Tval (Fields.Ventricle == 3)   +deg2rad(Fiber_info.AOTEPIRV) .*(1-Fields.Tval(Fields.Ventricle==3));   %OT interpolation
          
    %Endocardium orientation
          alpha_wall(Fields.Ventricle==-2)=deg2rad(Fiber_info.AENDOLV).*Fields.Tval_l(Fields.Ventricle==-2)  +deg2rad(Fiber_info.AOTENDOLV) .*(1-Fields.Tval_l(Fields.Ventricle==-2));  %OT interpolation
          alpha_wall(Fields.Ventricle==1 )=deg2rad(Fiber_info.AENDORV).*Fields.Tval(Fields.Ventricle==1)     +deg2rad(Fiber_info.AOTENDORV) .*(1-Fields.Tval(Fields.Ventricle==1));     %OT interpolation 
        
    %Myocardium orientation
   
        %RV
   
       indx=find(Fields.Ventricle==2);  %interpolacion RV
            alpha_wall_epi(indx)=deg2rad(Fiber_info.AEPIRV).*Fields.Tval (indx)   +deg2rad(Fiber_info.AOTEPIRV) .*(1-Fields.Tval(indx));  %OT interpolation
            alpha_wall_epi(indx)=alpha_wall_epi(indx).*(1-Fields.T_sept(indx))+Fields.T_sept(indx).*(Fiber_info.Septum_angleRV+alpha_wall_epi(indx)); %Septum interpolation
            alpha_wall_endo(indx)=deg2rad(Fiber_info.AENDORV).*Fields.Tval(indx)  +deg2rad(Fiber_info.AOTENDORV).*(1-Fields.Tval(indx));  %OT interpolation 
       
        alpha_wall(indx)=alpha_wall_endo(indx).*(1-d(indx))+alpha_wall_epi(indx).*(d(indx));   %linear transmural interpolation       
       
        %LV
   
       indx=find(Fields.Ventricle==-1); % interpolacion LV
            alpha_wall_epi(indx)=deg2rad(Fiber_info.AEPILV).*Fields.Tval_l (indx)   +deg2rad(Fiber_info.AOTEPILV) .*(1-Fields.Tval_l(indx));  %OT interpolation
         %interpolation of septum comes later         
            alpha_wall_endo(indx)=deg2rad(Fiber_info.AENDOLV).*Fields.Tval_l(indx)  +deg2rad(Fiber_info.AOTENDOLV).*(1-Fields.Tval_l(indx));  %OT interpolation
        
        alpha_wall(indx)=alpha_wall_endo(indx).*(1-d2(indx))+alpha_wall_epi(indx).*(d2(indx));   %linear transmural interpolation
          
%% 
        if Fiber_info.Interpolation_in_septum==1
    %    Correction of the septum points
    %     this step is neccesary in order to achieve a continuity in the septum.
    %     We can choose to modify the rv septal points or lv septal points.
    %     Here we modify the LV points (more points)
    %     For calculating this "angle of correction", we have to calculate
    %     individually for each septal point, and compare the e1 vector of the
    %     LV with the e1 vector of the RV 

    %(to improve)
            points_sept_all=find(Fields.T_sept>0.15);
            points_sept_RV=intersect(find(Fields.Ventricle==2),points_sept_all);
            points_sept_LV=intersect(find(Fields.Ventricle==-1 | Fields.Ventricle==-3),points_sept_all);

        if meshformat=="cut"
            indx=points_sept_RV; % interpolation RV
                alpha_wall_epi(indx)=deg2rad(Fiber_info.AEPIRV);  
                alpha_wall_epi(indx)=alpha_wall_epi(indx).*(1-Fields.T_sept(indx))+Fields.T_sept(indx).*(Fiber_info.Septum_angleRV); %Septum interpolation
                alpha_wall_endo(indx)=deg2rad(Fiber_info.AENDORV);  
                alpha_wall(indx)=alpha_wall_endo(indx).*(1-d(indx))+alpha_wall_epi(indx).*(d(indx)); 

            Septum_angleLV=Fiber_info.Discontinuity_angle;  

        else
            TR_sept=delaunayTriangulation(v(points_sept_RV,:));
            [N_sept,~]=nearestNeighbor(TR_sept,v(points_sept_LV,:));
            e1_rvsept=e1(points_sept_RV,:);
            e1_rvsept_ref=e1_rvsept(N_sept,:);
            e1_lvsept=e1(points_sept_LV,:);
            prod_vect=cross(e1_rvsept_ref,e1_lvsept,2);
            angle_correct=atan2(sqrt(sum(prod_vect.^2,2)),dot( e1_rvsept_ref,e1_lvsept,2));%angle between both e1 vectors

            Septum_angleLV=angle_correct+Fiber_info.Discontinuity_angle; 

        end
            
            indx=points_sept_LV; % interpolacion LV
                alpha_wall_epi(indx)=deg2rad(Fiber_info.AEPILV).*Fields.Tval_l (indx)   +deg2rad(Fiber_info.AOTEPILV) .*(1-Fields.Tval_l(indx));  %OT interpolation
                alpha_wall_epi(indx)=alpha_wall_epi(indx).*(1-Fields.T_sept(indx))+Fields.T_sept(indx).*(Septum_angleLV+Fiber_info.Septum_angleRV); %Septum interpolation
                alpha_wall_endo(indx)=deg2rad(Fiber_info.AENDOLV).*Fields.Tval_l(indx)  +deg2rad(Fiber_info.AOTENDOLV).*(1-Fields.Tval_l(indx));  %OT interpolation
                alpha_wall(indx)=alpha_wall_endo(indx).*(1-d2(indx))+alpha_wall_epi(indx).*(d2(indx)); 
        end
          
                        
             %% beta angle
                 
          beta_epiLVext =deg2rad(180-Fiber_info.beta_epiLV);
          beta_epiLVint =deg2rad(180+Fiber_info.beta_epiLV);
          beta_wall=Fiber_info.beta_endoRV*(1-d)+Fiber_info.beta_epiRV*(d);   %all heart (RV included)
          beta_wall(Fields.Tphi<=-1)=Fiber_info.beta_endoLV*(1-d(Fields.Tphi<=-1))+ beta_epiLVint*(d(Fields.Tphi<=-1));
          beta_wall(Fields.Tphi>-1 & Fields.Tphi<0 | Fields.Ventricle==-3 )=Fiber_info.beta_endoLV*(1-d(Fields.Tphi>-1 & Fields.Tphi<0 | Fields.Ventricle==-3 ))+beta_epiLVext*(d(Fields.Tphi>-1 & Fields.Tphi<0 | Fields.Ventricle==-3 ));
          
          beta_wall=beta_wall.*abs(Fields.T_sept-1)+deg2rad(180).*Fields.T_sept; %interpolacion septo
          beta_wall(Fields.Ventricle>0)=beta_wall(Fields.Ventricle>0).*Fields.Tval(Fields.Ventricle>0)+deg2rad(180).*(1-Fields.Tval(Fields.Ventricle>0));%interpolacion OT
          beta_wall(Fields.Ventricle<0)=beta_wall(Fields.Ventricle<0).*Fields.Tval_l(Fields.Ventricle<0)+deg2rad(180).*(1-Fields.Tval_l(Fields.Ventricle<0));%interpolacion OT
          
          beta_wall=beta_wall*Fiber_info.beta;  

          %% Rotation              
  %create matrix 3x3xN with vectors e0 e1 y e2 (column)
      Q=cat(3,e0,e1,e2); 
      Q=permute(Q,[3 2 1]);
     
    
   %rotation over e2
      axis=e2;
      theta=alpha_wall;
      Rot_axis=[cos(theta)+(axis(:,1).^2).*(1-cos(theta)) axis(:,1).*axis(:,2).*(1-cos(theta))-axis(:,3).*sin(theta) axis(:,1).*axis(:,3).*(1-cos(theta))+axis(:,2).*sin(theta);
             axis(:,2).*axis(:,1).*(1-cos(theta))+axis(:,3).*sin(theta) cos(theta)+(axis(:,2).^2).*(1-cos(theta)) axis(:,2).*axis(:,3).*(1-cos(theta))-axis(:,1).*sin(theta);
             axis(:,3).*axis(:,1).*(1-cos(theta))-axis(:,2).*sin(theta) axis(:,3).*axis(:,2).*(1-cos(theta))+axis(:,1).*sin(theta) cos(theta)+(axis(:,3).^2).*(1-cos(theta))];

        R=reshape(Rot_axis',3,N_points,3);
     R=permute(R,[3 1 2]);
     
      QX=zeros(size(R));
      for indx=1:length(alpha_wall)
          QX(:,:,indx)= mtimes(Q(:,:,indx), R(:,:,indx));
      end
  % Second rotation over ef
      axis2=QX(2,:,:);
      theta2=beta_wall;
      axis=reshape(axis2,3,N_points)';
      Rot_axis2=[cos(theta2)+(axis(:,1).^2).*(1-cos(theta2)) axis(:,1).*axis(:,2).*(1-cos(theta2))-axis(:,3).*sin(theta2) axis(:,1).*axis(:,3).*(1-cos(theta2))+axis(:,2).*sin(theta2);
             axis(:,2).*axis(:,1).*(1-cos(theta2))+axis(:,3).*sin(theta2) cos(theta2)+(axis(:,2).^2).*(1-cos(theta2)) axis(:,2).*axis(:,3).*(1-cos(theta2))-axis(:,1).*sin(theta2);
             axis(:,3).*axis(:,1).*(1-cos(theta2))-axis(:,2).*sin(theta2) axis(:,3).*axis(:,2).*(1-cos(theta2))+axis(:,1).*sin(theta2) cos(theta2)+(axis(:,3).^2).*(1-cos(theta2))];

        
    R2=reshape(Rot_axis2',3,N_points,3);
    R2=permute(R2,[3 1 2]);   
  
    QX2=zeros(size(R));

    for indx=1:length(alpha_wall)
        QX2(:,:,indx)= mtimes(QX(:,:,indx), R2(:,:,indx));
    end

    FX=QX2(1,:,:);
    Fields.F=reshape(FX,3,N_points)';

    FY=QX2(2,:,:);
    Fields.F_N=reshape(FY,3,N_points)';

    FZ=QX2(3,:,:);
    Fields.F_S=reshape(FZ,3,N_points)';

% clear -regexp ^Q;
clear -regexp ^R;                   
        %% POSTPROCESSING (valves, apex presentation...)  
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Post-processing');


% transmurality

Fields.d3=abs(abs(Fields.Tphi)-1); 
Fields.d3(Fields.Tphi<0)=abs(abs(Fields.Tphi(Fields.Tphi<0))/2-1);
Fields.d3=(0-1)*(Fields.d3-min(Fields.d3))./(max(Fields.d3)-min(Fields.d3))+(1);


 %%               

disp('Coordinate calculation');

%apex2base
%Geodesic distance calculation
%Distance calculation--> heat method
FieldsC.Tpsi_bi(FieldsC.Tpsi_bi>0.99)=1.0;%
Apex_points_id=find(FieldsC.Tpsi_bi==1.0);
[phi] = heat_geodesic(v2,double(f2(:,2:end)+1),double(Apex_points_id),[],'BoundaryConditions','robin','Legacy',1);

%reescale distances between 1 and 2
FieldsC.apex2base=(phi-min(phi))./(max(phi)-min(phi));  %just geodesic distance
%a different approuch could be projecting the distances in the Z axis

FieldsC.apex_2_base=(1-2)*(FieldsC.apex2base-min(FieldsC.apex2base))./(max(FieldsC.apex2base)-min(FieldsC.apex2base))+2;


%% Material (Tissue conductivity) especification

if meshformat=="cut"
   Fields.Plug_points=zeros(size(Fields.Tpsi_bi));
   Fields.Plug_points((Fields.Tpsi_bi <=0.0005 & Fields.Ventricle>0) | (Fields.Tpsi_bi <=0.0005 & Fields.Ventricle<0))=1;
   Fields.Plug_points(abs(Fields.Tphi)<eps)=0;
else
   Fields.Plug_points=zeros(size(Fields.TpsiR));
   Fields.Plug_points((Fields.TpsiR<=0.0005 & Fields.Ventricle>0) | (Fields.TpsiL<=0.0005 & Fields.Ventricle<0))=1;
   Fields.Plug_points(abs(Fields.Tphi)<eps)=0;
end

  Plug_points_id=find(Fields.Plug_points); 
  [Plug_tetra0]=ismember(f(:,2:end)+1,Plug_points_id);
  Fields.Plug_tetra=max(Plug_tetra0*1,[],2)+1;

%% 
                 
%Epiendo transmurality
Fields.Epiendo=ones(size(Fields.d3)).*2;
Fields.Epiendo(Fields.d3>=(1-epiendo(1)./100))=1;
Fields.Epiendo(Fields.d3<(epiendo(3)./100))=3;


%RV endocardium as epi
Fields.Epiendo3=ones(size(Fields.Tphi3)).*2;
Fields.Epiendo3(Fields.Tphi3>=(1-epiendoRV(1)./100))=1;
Fields.Epiendo3(Fields.Tphi3<(epiendoRV(3)./100))=3;


% %%lvrv
% body=Fields.Tphi(elem)./abs(Fields.Tphi(elem));
% body1=mode(body,2);
% 
%   %in order to classify NaN, we look for the neighbors
%  TR_mesh=triangulation(double(elem),v2);
%   for iter=1:10
%       zeros_body=find(isnan(body1));
%       Neighbors_zeros = neighbors(TR_mesh,zeros_body);
%       Neighbors_val=body1(mode(Neighbors_zeros,2));
%       body1(zeros_body)=Neighbors_val;
%   end

  
lvrv=zeros(size(label_tetra));
lvrv(body1>0)=1;
FieldsC.lvrv=lvrv;

%calculate ventricular gradients 
[FieldsC.a2b,FieldsC.r2l,FieldsC.a2p,a2b_vector,r2l_vector,a2p_vector]=ventricular_directions(v2,pto,car,Fid2,face2,labelf2);
  %write_vtk_rbm(strcat('Grads.vtk'),length(v2),v2,length(f2),f2,FieldsC.a2b,FieldsC.r2l,FieldsC.a2p); %write mesh 

%cobiveco coord
[FieldsC.r,FieldsC.a2b_cobi,FieldsC.tm_cobi,FieldsC.lvrv_cobi,apex_cobi_id]=cobi_coords200_UKBB(node,labelp,FieldsC.a2b,FieldsC.a2p,FieldsC.T_sept_cobi,elem,label,labelf3,face,face3,Fid,FieldsC.Ventricle,FieldsC.Tv_cobi,L,GT);
[FieldsC.a2b_uvc]=heat_method3(node,elem,apex_cobi_id,FieldsC.DTpsi_bi);
FieldsC.rsin=sin(2*pi*FieldsC.r);
FieldsC.rcos=cos(2*pi*FieldsC.r);

FieldsC.r2l_geo=heat_geodesic(v2,double(f2(:,2:end)+1),find(FieldsC.T_sept_uvc>0.95) ,[],'BoundaryConditions','robin','Legacy',1);

    r2l_geo(FieldsC.Tv_cobi<=0.5)= normalize(FieldsC.r2l_geo(FieldsC.Tv_cobi<=0.5), 'range', [0 1]);
    r2l_geo(FieldsC.Tv_cobi>0.5)=normalize(FieldsC.r2l_geo(FieldsC.Tv_cobi>0.5)*-1, 'range', [-1 0]);


%cut geometry a2b equivalence
if meshformat=="cut"
    FieldsC.a2b_cut=FieldsC.a2b_uvc;
else
    indValve=(labelp==9 | labelp==10) |(labelp==13 | labelp==14);
    Threshold_val=min(FieldsC.a2b(indValve));
    FieldsC.a2b_cut=FieldsC.a2b_uvc;
    FieldsC.a2b_cut=(FieldsC.a2b_cut-min(FieldsC.a2b_cut))./(Threshold_val-min(FieldsC.a2b_cut));
end

%pericardium
disp('Pericardium generation');
centroid_surf=meshcentroid(v2,double(face));
warning('off','all')
TRinv=triangulation(double(face),v2);
warning('on','all')
[NN_surf,~]=nearestNeighbor(TRinv,centroid_surf);
label_A2B=FieldsC.a2b_uvc(NN_surf);

label_set=label_closed;
label_set(label_closed==3 | label_closed==6 | label_closed==18 | label_closed==9 | label_closed==14)=2;
label_set(label_closed==2  | label_closed==10 | label_closed==13 )=3;
label_set( label_closed==12 |label_set>18)=1;

FieldsC.label_set=label_set;
label_fine=label;
label_set2=label_set;
label_set2(label_set==1)=4;
label_set2((label_A2B <pericardiumlevel & label_set==1) )=1;
FieldsC.label_set2=label_set2;
%write_vtk_surf('labelsP.vtk',v2,face,FieldsC.label_set2);

cd (directoryResults)

%% interpolation to finer mesh
disp('final interpolation');
 if largenumber==1
     %vertex
     apex_2_basepointsTetra=FieldsC.apex_2_base(points_Tetra);
     FieldsF.apex_2_base=dot(bar,apex_2_basepointsTetra,2);

     a2bpointsTetra=FieldsC.a2b(points_Tetra);
     FieldsF.a2b=dot(bar,a2bpointsTetra,2);

     r2lpointsTetra=FieldsC.r2l(points_Tetra);
     FieldsF.r2l=dot(bar,r2lpointsTetra,2);

     a2ppointsTetra=FieldsC.a2p(points_Tetra);
     FieldsF.a2p=dot(bar,a2ppointsTetra,2);

     rsinpointsTetra=FieldsC.rsin(points_Tetra);
     FieldsF.rsin=dot(bar,rsinpointsTetra,2);

     rcospointsTetra=FieldsC.rcos(points_Tetra);
     FieldsF.rcos=dot(bar,rcospointsTetra,2);

     r2l_geopointsTetra=FieldsC.r2l_geo(points_Tetra);
     FieldsF.r2l_geo=dot(bar,r2l_geopointsTetra,2);

     a2b_cobipointsTetra=FieldsC.a2b_cobi(points_Tetra);
     FieldsF.a2b_cobi=dot(bar,a2b_cobipointsTetra,2);

     tm_cobipointsTetra=FieldsC.tm_cobi(points_Tetra);
     FieldsF.tm_cobi=dot(bar,tm_cobipointsTetra,2);

     a2b_uvcpointsTetra=FieldsC.a2b_uvc(points_Tetra);
     FieldsF.a2b_uvc=dot(bar,a2b_uvcpointsTetra,2);

     a2b_cutpointsTetra=FieldsC.a2b_cut(points_Tetra);
     FieldsF.a2b_cut=dot(bar,a2b_cutpointsTetra,2);

        %projected cobiveco

     %% 

     %faces

     centroids_fine=meshcentroid(node_vox,face_fine);
     TR_coarse=delaunayTriangulation(centroid);
     [NN_faces,~]=nearestNeighbor(TR_coarse,centroids_fine);
     FieldsF.label_set2=FieldsC.label_set2(NN_faces);
     FieldsF.label_fine=label_closed(NN_faces);

     %elements
     centroids_fine_elem=meshcentroid(node_vox,f(:,2:end)+1);
     TR_elem=delaunayTriangulation(centroid_tetra);
     [NN_elem,~]=nearestNeighbor(TR_elem,centroids_fine_elem);
     FieldsF.lvrv=lvrv(NN_elem);

     FieldsF.lvrv_cobi=ones(size(FieldsF.Ventricle));
     FieldsF.lvrv_cobi(FieldsF.Tv_cobi<=0.5)=0;

     %Transform r coordinate

      FieldsF.r=atan2( FieldsF.rsin, FieldsF.rcos)./(2*pi);
      FieldsF.r( FieldsF.rsin<0)=atan2( FieldsF.rsin( FieldsF.rsin<0), FieldsF.rcos( FieldsF.rsin<0))./(2*pi)+1;
      FieldsF.r(isnan( FieldsF.r))=-1;
      FieldsF.a2b_cut( FieldsF.r==-1)=-1;
      FieldsF.a2b_cobi(isnan( FieldsF.a2b_cobi))=-1;
      %% create ensi cases
      disp('Saving data');
        %%aha
     [ FieldsF.aha]= aha_segments_v3(v,f, FieldsF.Ventricle, FieldsF.r,pto,car+1,Fid);

    mkdir(fullfile('ensi_Fine_',num2str(case_number)))
    directoryFine=fullfile(directoryResults,'ensi_Fine_',num2str(case_number));
    cd (directoryFine)
    save_ensi_UKBB(v,f, FieldsF.Ventricle, Fields.d3,FieldsF.Tphi3,FieldsF.tm_cobi,Fields.Epiendo,Fields.Epiendo3,FieldsF.a2b_uvc,FieldsF.a2b_cobi,FieldsF.r,FieldsF.lvrv_cobi,FieldsF.r2l_geo,FieldsF.a2b,FieldsF.r2l,FieldsF.a2p,Fields.F,Fields.F_S,Fields.F_N,FieldsF.apex_2_base,FieldsF.aha,Fields.Plug_tetra);
    Tphi3=FieldsF.Tphi3; d3=Fields.Tphi3;Plug_points=Fields.Plug_points;Tphi=FieldsF.Tphi;Tphi_bi=FieldsF.Tphi_bi;F=Fields.F;F_S=Fields.F_S;F_N=Fields.F_N;
    Epiendo=Fields.Epiendo; Epiendo3=Fields.Epiendo3; Ventricle=FieldsF.Ventricle; Plug_tetra=Fields.Plug_tetra;apex_2_base=FieldsF.apex_2_base; label_set2=FieldsF.label_set2;
    r=FieldsF.r; aha=FieldsF.aha;tm_cobi=FieldsF.tm_cobi;a2b_uvc=FieldsF.a2b_uvc;r2l_geo=FieldsF.r2l_geo;lvrv_cobi=FieldsF.lvrv_cobi;a2b_cobi=FieldsF.a2b_cobi;a2b=FieldsF.a2b;a2p=FieldsF.a2p;r2l=FieldsF.r2l; a2b_cut=FieldsF.a2b_cut;lvrv=FieldsF.lvrv;label_fine=FieldsF.label_fine;
    save('Case_Fine','v','f','Tphi3','d3','Plug_points','Tphi','Tphi_bi','Epiendo','Epiendo3','Ventricle','F','F_N','F_S','Plug_tetra','label','label_set','label_set2','apex_2_base','lvrv','NN_surf','r','tm_cobi','a2b_uvc','a2b_cobi','lvrv_cobi','r2l_geo','a2b','a2p','r2l','a2b_cut','aha','a2b_vector','r2l_vector','a2p_vector','label_fine');
    CSVFilesgeneration_f_UKBB([num2str(case_number),'_',name(1:end-4)],'Case_fine.mat');

     cd ..
     %
     %write_vtk_rbm(strcat('Info.vtk'),length(v),v,length(f),f,Ventricle,Epiendo,Tphi3,F); %write mesh with fibers Long

     %fine to coarse mesh IDs
     TR_fine=triangulation(double(f(:,2:end)+1),v);
     [NN,~]=nearestNeighbor(TR_fine,v2);
     %fine to coarse elems
     TR_elem_f=delaunayTriangulation(centroids_fine_elem);
     [NN_elem_c,~]=nearestNeighbor(TR_elem_f,centroid_tetra);
     %coarse to fine mesh IDs
     TR_coarse=triangulation(double(f2(:,2:end)+1),v2);
     [NN_c2f,~]=nearestNeighbor(TR_coarse,v);
    mkdir(strcat('ensi',num2str(case_number)))
    directoryFine=fullfile(directoryResults,'ensi',num2str(case_number));
    cd (directoryFine)

      %interp
     Plug_points=Plug_points(NN); Plug_points_id=find(Plug_points);% [Plug_tetra0]=ismember(f2(:,2:end)+1,Plug_points_id); Plug_tetra=max(Plug_tetra0*1,[],2)+1; 
     F=F(NN,:); F_S=F_S(NN,:); F_N=F_N(NN,:);aha=FieldsF.aha(NN,:);
     Plug_tetra=Fields.Plug_tetra(NN_elem_c,:);
    % transmurality Coarse
        FieldsC.d3=abs(abs(FieldsC.Tphi)-1); 
        FieldsC.d3(FieldsC.Tphi<0)=abs(abs(FieldsC.Tphi(FieldsC.Tphi<0))/2-1);
        FieldsC.d3=(0-1)*(FieldsC.d3-min(FieldsC.d3))./(max(FieldsC.d3)-min(FieldsC.d3))+(1);
        %Epiendo transmurality Coarse
        FieldsC.Epiendo=ones(size(FieldsC.d3)).*2;
        FieldsC.Epiendo(FieldsC.d3>=(1-epiendo(1)./100))=1;
        FieldsC.Epiendo(FieldsC.d3<(epiendo(3)./100))=3;   
        %RV endocardium as epi Coarse
        FieldsC.Epiendo3=ones(size(FieldsC.Tphi3)).*2;
        FieldsC.Epiendo3(FieldsC.Tphi3>=(1-epiendoRV(1)./100))=1;
        FieldsC.Epiendo3(FieldsC.Tphi3<(epiendoRV(3)./100))=3;
      save_ensi_UKBB(v2,f2,FieldsC.Ventricle,FieldsC.d3,FieldsC.Tphi3,FieldsC.tm_cobi,FieldsC.Epiendo,FieldsC.Epiendo3,FieldsC.a2b_uvc,FieldsC.a2b_cobi,FieldsC.r,FieldsC.lvrv_cobi,FieldsC.r2l_geo,FieldsC.a2b,FieldsC.r2l,FieldsC.a2p,F,F_S,F_N,FieldsC.apex_2_base,aha,Plug_tetra);

     save('IDs_rel','NN','NN_c2f');
     %coarsening:
     v=v2;
     f=f2;
     label_fine=label;
     lvrv=FieldsC.lvrv;
      FieldsC.r(isnan( FieldsC.r))=-1;
      FieldsC.a2b_cut( FieldsC.r==-1)=-1;
      FieldsC.a2b_cobi(isnan( FieldsC.a2b_cobi))=-1;

    Tphi3=FieldsC.Tphi3; d3=FieldsC.Tphi3;Tphi=FieldsC.Tphi;Tphi_bi=FieldsC.Tphi_bi;
    Epiendo=FieldsC.Epiendo; Epiendo3=FieldsC.Epiendo3; Ventricle=FieldsC.Ventricle; apex_2_base=FieldsC.apex_2_base; label_set2=FieldsC.label_set2;
    r=FieldsC.r;tm_cobi=FieldsC.tm_cobi;a2b_uvc=FieldsC.a2b_uvc;r2l_geo=FieldsC.r2l_geo;lvrv_cobi=FieldsC.lvrv_cobi;a2b_cobi=FieldsC.a2b_cobi;a2b=FieldsC.a2b;a2p=FieldsC.a2p;r2l=FieldsC.r2l;a2b_cut=FieldsC.a2b_cut;

     save('Case_coarse','v','f','Tphi3','d3','Plug_points','Tphi','Tphi_bi','Epiendo','Epiendo3','Ventricle','F','F_N','F_S','Plug_tetra','label','label_set','label_set2','apex_2_base','lvrv','NN_surf','r','tm_cobi','a2b_uvc','a2b_cobi','lvrv_cobi','r2l_geo','a2b','a2p','r2l','a2b_cut','aha','a2b_vector','r2l_vector','a2p_vector','label_fine');
     CSVFilesgeneration_f_UKBB(strcat(num2str(case_number),'_coarse'),'Case_coarse.mat');

 else  
disp('Saving data');
     FieldsC.r(isnan(FieldsC.r))=-1;
     FieldsC.a2b_cut(FieldsC.r==-1)=-1;
     FieldsC.a2b_cobi(isnan( FieldsC.a2b_cobi))=-1;

     %%aha
     [Fields.aha]= aha_segments_v3(v,f,FieldsC.Ventricle,FieldsC.r,pto,car+1,Fid);
     mkdir('ensi')
     directoryFine=fullfile(directoryResults,'ensi');
     cd (directoryFine)
     save_ensi_UKBB(v,f,Fields.Ventricle,Fields.d3,Fields.Tphi3,FieldsC.tm_cobi,Fields.Epiendo,Fields.Epiendo3,FieldsC.a2b_uvc,FieldsC.a2b_cobi,FieldsC.r,FieldsC.lvrv_cobi,FieldsC.r2l_geo,FieldsC.a2b,FieldsC.r2l,FieldsC.a2p,Fields.F,Fields.F_S,Fields.F_N,FieldsC.apex_2_base,Fields.aha,Fields.Plug_tetra);
%% save data

    Tphi3=Fields.Tphi3; d3=Fields.Tphi3;Plug_points=Fields.Plug_points;Tphi=Fields.Tphi;Tphi_bi=Fields.Tphi_bi;F=Fields.F;F_S=Fields.F_S;F_N=Fields.F_N;
    Epiendo=Fields.Epiendo; Epiendo3=Fields.Epiendo3; Ventricle=Fields.Ventricle; Plug_tetra=Fields.Plug_tetra;apex_2_base=FieldsC.apex_2_base; label_set2=FieldsC.label_set2;
    r=FieldsC.r; aha=Fields.aha;tm_cobi=FieldsC.tm_cobi;a2b_uvc=FieldsC.a2b_uvc;r2l_geo=FieldsC.r2l_geo;lvrv_cobi=FieldsC.lvrv_cobi;a2b_cobi=FieldsC.a2b_cobi;a2b=FieldsC.a2b;a2p=FieldsC.a2p;r2l=FieldsC.r2l;a2b_cut=FieldsC.a2b_cut;

    save('Case','v','f','Tphi3','d3','Plug_points','Tphi','Tphi_bi','Epiendo','Epiendo3','Ventricle','F','F_N','F_S','Plug_tetra','label','label_set','label_set2','apex_2_base','lvrv','NN_surf','r','tm_cobi','a2b_uvc','a2b_cobi','lvrv_cobi','r2l_geo','a2b','a2p','r2l','aha','a2b_vector','r2l_vector','a2p_vector','label_fine','a2b_cut');
    CSVFilesgeneration_f_UKBB(name(1:end-4),'Case.mat');
    %
    %write_vtk_rbm(strcat('Infoc.vtk'),length(v),v,length(f),f,Ventricle,Epiendo,Tphi3,F); %write mesh with fibers Long

 end

cd (directory)

%% functions
    function pointGrad=elemGrad2pointGrad(elemGrad,TR)
             %project the gradients calculated in the elements to the vertexes
                tri_id= vertexAttachments(TR,(1:length(TR.Points))');
                pointGrad=zeros(length(tri_id),3);
                for i=1:length(tri_id)            
                    pointGrad(i,:)=mean(elemGrad(tri_id{i,1},:),1);            
                end
    
    end
    function  Sorted_Con=sortconnectivities(connectivities_orig)
            %reorder mesh connectivities
                [a1, ~, a2] = unique(connectivities_orig);
                Anew=1:length(a1);
                Surfaces0 = Anew(a2);
                Sorted_Con= reshape(Surfaces0, size(connectivities_orig));
    end

    function values_interp = baryInterpIgnoreNaN(bar, vertex_values)
    % BARYINTERPIGNORENAN - Barycentric interpolation that ignores NaN vertex values
    %
    % Inputs:
    %   bar           - [N x 4] barycentric coordinates
    %   vertex_values - [N x 4] values at the tetrahedron vertices (may include NaNs)
    %
    % Output:
    %   values_interp - [N x 1] interpolated values, ignoring NaNs
    
        valid = ~isnan(vertex_values);                    % Logical mask of valid entries
        vertex_values(~valid) = 0;                        % Set NaNs to zero
        weighted_bar = bar .* valid;                      % Zero out weights for NaNs
        norm_factor = sum(weighted_bar, 2);               % Sum of valid weights
        norm_bar = weighted_bar ./ norm_factor;           % Renormalize weights
        values_interp = sum(norm_bar .* vertex_values, 2);% Interpolated result
    end

end

