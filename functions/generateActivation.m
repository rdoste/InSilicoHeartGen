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


function [point_time_new,IDs_endo]=generateActivation(name,case_final_name,reference_folder,Folder_name,tag,BCL)
    
    %generateActivation('Case_Fine.mat',"generic")
    %tag--> "personalised" or "generic"
    %if is personalised, reads csv file from specified folder
    %if is generic, reads projects activation from generic mesh (CR05)

    % tag="personalised";
    % case_final_name="Case_Fine.mat";
    % BCL=0.8;
    % 
    if tag=="generic"
    
        %read files
        case_final=case_final_name;


        case_orig0=load(fullfile(reference_folder,'Reference.mat'));
        case_orig=case_orig0.Reference;
        namefile=fullfile(reference_folder,'Reference.stimuli');
        fid = fopen(namefile);
        point_data=textscan(fid, '%d %f %f %f %f %f ','HeaderLines',1); 
        fclose(fid);
        point_ID=double(point_data{1,1}(:));
        current=point_data{1,2}(:);
        point_time=point_data{1,3}(:);
        duration=point_data{1,4}(:);
        beats=point_data{1,5}(:);
        bcl=point_data{1,6}(:);
        write_vtk_points('Activation_points_orig.vtk',case_orig.v,point_ID,point_time);

         %%%%%%%%%manual correction for cut geometries%%%%%%%%%%%%%%%%%%%
        %lower basal plane in the final mesh
        case_final.a2b_cobi(case_final.a2b_cobi==-1)=NaN;
        case_orig.a2b_cobi( case_orig.a2b_cobi==-1)=NaN;
        cutoff_value=1;   % 1 normal
  
         %%%%%%%
        
        
        %look for nearest points
        endopoints_orig_LV=intersect(find(case_orig.Ventricle==-2 & case_orig.a2b_cobi<=cutoff_value),point_ID);
        endopoints_final_LV=find(case_final.Ventricle==-2 & case_final.a2b_cobi<=cutoff_value);
        
        
         [k_L,distL] = dsearchn([case_orig.a2b_cobi(endopoints_orig_LV),case_orig.r(endopoints_orig_LV)],[case_final.a2b_cobi(endopoints_final_LV),case_final.r(endopoints_final_LV)]);
        
        endopoints_orig_RV =intersect(find(case_orig.Ventricle==1 & case_orig.a2b_cobi<=cutoff_value),point_ID);
        endopoints_final_RV=find(case_final.Ventricle==1 & case_final.a2b_cobi<=cutoff_value );
        
         [k_R,distR] = dsearchn([case_orig.a2b_cobi(endopoints_orig_RV),case_orig.r(endopoints_orig_RV)],[case_final.a2b_cobi(endopoints_final_RV),case_final.r(endopoints_final_RV)]);
        
        
        % Activation projection
        
        point_time_full=zeros(length(case_orig.v),1);
        point_time_full(point_ID)=point_time;
        point_time_LV=point_time_full(endopoints_orig_LV);
        point_time_new1=point_time_LV(k_L);
        point_time_RV=point_time_full(endopoints_orig_RV);
        point_time_new2=point_time_RV(k_R);


        %% map the basal points using gradients
        endopoints_orig_basal_LV=intersect(find(case_orig.Ventricle==-2 & isnan(case_orig.a2b_cobi)),point_ID);
        endopoints_final_basal_LV=intersect(find(case_final.Ventricle==-2 & isnan(case_final.a2b_cobi)),find(case_final.Plug_points==0));
        [k_B_LV,distb] = dsearchn([case_orig.a2b_uvc(endopoints_orig_basal_LV),case_orig.r2l_geo(endopoints_orig_basal_LV),case_orig.a2p(endopoints_orig_basal_LV)],[case_final.a2b_uvc(endopoints_final_basal_LV),case_final.r2l_geo(endopoints_final_basal_LV),case_final.a2p(endopoints_final_basal_LV)]);

        endopoints_orig_basal_RV=intersect(find(case_orig.Ventricle==1 & isnan(case_orig.a2b_cobi)),point_ID);
        endopoints_final_basal_RV=intersect(find(case_final.Ventricle==1 & isnan(case_final.a2b_cobi)),find(case_final.Plug_points==0));
        [k_B_RV,distbR] = dsearchn([case_orig.a2b_uvc(endopoints_orig_basal_RV),case_orig.r2l_geo(endopoints_orig_basal_RV),case_orig.a2p(endopoints_orig_basal_RV)],[case_final.a2b_uvc(endopoints_final_basal_RV),case_final.r2l_geo(endopoints_final_basal_RV),case_final.a2p(endopoints_final_basal_RV)]);

        point_time_LV_basal=point_time_full(endopoints_orig_basal_LV);
        point_time_new1_basal=point_time_LV_basal(k_B_LV);
        point_time_RV_basal=point_time_full(endopoints_orig_basal_RV);
        point_time_new2_basal=point_time_RV_basal(k_B_RV);

%% 
        point_time_new=[point_time_new1; point_time_new2;point_time_new1_basal;point_time_new2_basal] ;
        point_ID_new=[endopoints_final_LV;endopoints_final_RV;endopoints_final_basal_LV;endopoints_final_basal_RV];
        
        Act_modif=[point_ID_new,repmat(-200,length(point_ID_new),1),point_time_new,repmat(0.001,length(point_ID_new),1),repmat(BCL*10,length(point_ID_new),1),repmat(BCL,length(point_ID_new),1)];
        
        
        %Sort in ascinding order
        [~,idu] = unique(Act_modif(:,1));
        A= Act_modif(idu,:);
        
        cd (Folder_name)
        write_vtk_points('Activation_points_final.vtk',case_final.v,point_ID_new,point_time_new);         
        
        fid = fopen(strcat(name,'.stimuli'),'w');
        fprintf(fid,'%d %d %.9f %.3f %d %f \n',A');  
        fclose(fid);
        IDs_endo=point_ID_new;
        point_time_new=point_time_new;
    
    
    elseif tag=="personalised"

        %read files
        case_final=load(case_final_name);
        case_corse=load('Case.mat');
        IDs_correspondance=load("IDs_rel.mat");
        cd ('InferenceResults')
        LATFile=dir('*lat.csv');
        LAT_coarse=readmatrix(LATFile.name);
        cd ..

        LAT_fine=LAT_coarse(IDs_correspondance.NN_c2f);
        IDs_endo=find(case_final.Ventricle==1 | case_final.Ventricle==-2);
        LAT_fine_endo=LAT_fine(IDs_endo)./1000; %Endo points and conversion to seconds
        Act=[IDs_endo,repmat(-200,length(IDs_endo),1),LAT_fine_endo,repmat(0.001,length(IDs_endo),1),repmat(8,length(LAT_fine_endo),1),repmat(0.8,length(LAT_fine_endo),1)];

        write_vtk_points('Activation_points_final_person.vtk',case_final.v,IDs_endo,LAT_fine_endo);
        
        
        fid = fopen('Final_person.stimuli','w');
        fprintf(fid,'%d %d %.9f %.3f %d %f \n',Act');
    
    else 
        error('wrong tag input')
    end
% end