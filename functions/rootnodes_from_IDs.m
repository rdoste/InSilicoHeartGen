function [rootnodes]=rootnodes_from_IDs(name_alya,case_final,reference_folder,AlyaFolder,nroots,BCL)
    

    
        %transform to center in monoalg
        case_final.v=[case_final.v(:,1)-min(case_final.v(:,1)),case_final.v(:,2)-min(case_final.v(:,2)),case_final.v(:,3)-min(case_final.v(:,3))];
       


        case_orig0=load(strcat(reference_folder,'\Reference.mat'));
        case_orig=case_orig0.Reference;
        original_roots=load(strcat(reference_folder,'\roots_IDs.mat'));
        point_ID=original_roots.roots_IDs(:,1);
        point_time=original_roots.roots_IDs(:,2);
         %%%%%%%%%manual correction for DTI004 projection%%%%%%%%%%%%%%%%%%%
        %lower basal plane in the final mesh
        case_final.a2b_cut(case_final.a2b_cut==-1)=NaN;
        case_orig.a2b_cut( case_orig.a2b_cut==-1)=NaN;
        cutoff_value=1;   % 1 normal
  
         %%%%%%%
        roots=[];
        
        %look for nearest points
        endopoints_orig_LV=intersect(find(case_orig.Ventricle==-2 & case_orig.a2b_cut<=cutoff_value),point_ID);
        endopoints_final_LV=find(case_final.Ventricle==-2 & case_final.a2b_cut<=cutoff_value);
        
        if ~isempty (endopoints_orig_LV)
         [k_L,distL] = dsearchn([case_final.a2b_cobi(endopoints_final_LV),case_final.r(endopoints_final_LV)],[case_orig.a2b_cobi(endopoints_orig_LV),case_orig.r(endopoints_orig_LV)]);
        
         [~,roots_original_position]=ismember(endopoints_orig_LV,point_ID);
         roots(roots_original_position)=endopoints_final_LV(k_L);
        end



        endopoints_orig_RV =intersect(find(case_orig.Ventricle==1 & case_orig.a2b_cut<=cutoff_value),point_ID);
        endopoints_final_RV=find(case_final.Ventricle==1 & case_final.a2b_cut<=cutoff_value );
        if ~isempty (endopoints_orig_RV)
         [k_R,distR] = dsearchn([case_final.a2b_cobi(endopoints_final_RV),case_final.r(endopoints_final_RV)],[case_orig.a2b_cobi(endopoints_orig_RV),case_orig.r(endopoints_orig_RV)]);
         [~,roots_original_position]=ismember(endopoints_orig_RV,point_ID);
         roots(roots_original_position)=endopoints_final_RV(k_R);        
        end
        
        % Activation projection
        



        %% map the basal points using gradients
        endopoints_orig_basal_LV=intersect(find(case_orig.Ventricle==-2 & isnan(case_orig.a2b_cut)),point_ID);
        endopoints_final_basal_LV=intersect(find(case_final.Ventricle==-2 & isnan(case_final.a2b_cut)),find(case_final.Plug_points==0));
        if ~isempty (endopoints_orig_basal_LV)
            [k_B_LV,distb] = dsearchn([case_final.a2b_uvc(endopoints_final_basal_LV),case_final.r2l_geo(endopoints_final_basal_LV),case_final.a2p(endopoints_final_basal_LV)],[case_orig.a2b_uvc(endopoints_orig_basal_LV),case_orig.r2l_geo(endopoints_orig_basal_LV),case_orig.a2p(endopoints_orig_basal_LV)]);
            [~,roots_original_position]=ismember(endopoints_orig_basal_LV,point_ID);
            roots(roots_original_position)=endopoints_final_basal_LV(k_B_LV);  
        end
        endopoints_orig_basal_RV=intersect(find(case_orig.Ventricle==1 & isnan(case_orig.a2b_cut)),point_ID);
        endopoints_final_basal_RV=intersect(find(case_final.Ventricle==1 & isnan(case_final.a2b_cut)),find(case_final.Plug_points==0));
        if ~isempty (endopoints_orig_basal_RV)
            [k_B_RV,distbR] = dsearchn([case_final.a2b_uvc(endopoints_final_basal_RV),case_final.r2l_geo(endopoints_final_basal_RV),case_final.a2p(endopoints_final_basal_RV)],[case_orig.a2b_uvc(endopoints_orig_basal_RV),case_orig.r2l_geo(endopoints_orig_basal_RV),case_orig.a2p(endopoints_orig_basal_RV)]);
            [~,roots_original_position]=ismember(endopoints_orig_basal_RV,point_ID);
            roots(roots_original_position)=endopoints_final_basal_RV(k_B_RV);  
        end



      rootnodes.coords=case_final.v(roots,:);
      rootnodes.IDs=roots;
      rootnodes.time=point_time;


      write_vtk_points('Activation_roots.vtk',case_final.v(roots,:),(1:nroots),point_time);
 
    
 end