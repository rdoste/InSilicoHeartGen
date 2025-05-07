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


function HexaFieldsGeneration_function_cells_v2(monodir,ALG,MeshCoarse,MeshHex,Data,tet_ID,bar,epiendo,epiendoRV,case_name)
%add original projected cobiveco coodinates

                    
   %% Interpolation for cube centers

 disp('Interpolation to hexahedral mesh');
      
                           
    node_vox=double(MeshHex.points);                 
    elem_vox=MeshHex.cells;
    elem=[ones(length(elem_vox),1)*4,elem_vox-1];  %remove in the future

    cube_centroid=ALG(:,1:3)/1000;

  
    % Fields interpolation
    
    points_Tetra=MeshCoarse.cells(tet_ID,:);
    
    % Interpolation
   
    TphipointsTetra=Data.Tphi(points_Tetra);
    Tphi=dot(bar,TphipointsTetra,2);
    
    TphibipointsTetra=Data.Tphi_bi(points_Tetra);
    Tphi_bi=dot(bar,TphibipointsTetra,2);

    TphicobipointsTetra=Data.tm_cobi(points_Tetra);
    tm_cobi=dot(bar,TphicobipointsTetra,2);


    Tphi3pointsTetra=Data.Tphi3(points_Tetra);
    Tphi3=dot(bar,Tphi3pointsTetra,2);
  
    Material=Data.Plug_points+1; %to set 1 as healthy tissue
    MaterialpointsTetra=Material(points_Tetra);
    Material=dot(bar,MaterialpointsTetra,2);



%Ventricle Tag interpolation
    
      Ventricle2=int8(Data.Ventricle(points_Tetra));
      Ventricle_raw=round(dot(bar,double(Ventricle2),2));%know lV vs RV 
      %rounded value
      Epiendo3=Ventricle2;%auxiliar
      %positive values
      Epiendo3(Ventricle2==-2)=1;
      Epiendo3(Ventricle2==-1)=2;
      Epiendo3(Ventricle2==-3)=3;
      %Ventricle values calculation
      Ventricle=round(dot(bar,double(Epiendo3),2));
      Ventricle(Ventricle_raw<0)=Ventricle(Ventricle_raw<0).*-1;            
      Epiendo3=int8(Ventricle);%auxiliar
      Ventricle(Epiendo3==-1)=-2;
      Ventricle(Epiendo3==-2)=-1;
clear -regexp pointsTetra$;                        
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Post-processing');

% transmurality

d3=abs(abs(Tphi)-1); 
d3(Tphi<0)=abs(abs(Tphi(Tphi<0))/2-1);


%% Material (Tissue conductivity) especification

    Material=Data.Plug_points+1; %to set 1 as healthy tissue
    MaterialpointsTetra=Material(points_Tetra);
    Material=max(MaterialpointsTetra,[],2);

  

%% 

%Epiendo transmurality
Epiendo=ones(size(d3)).*2;
Epiendo(d3>=(1-epiendo(1)./100))=1;
Epiendo(d3<(epiendo(3)./100))=3;


%RV endocardium as epi
Epiendo3=ones(size(Tphi3)).*2;
Epiendo3(Tphi3>=(1-epiendoRV(1)./100))=1;
Epiendo3(Tphi3<(epiendoRV(3)./100))=3;



%% interpolation to finer mesh
Data.r(Data.r==-1)=NaN;
rsin=sin(2*pi*Data.r);
rcos=cos(2*pi*Data.r);

     %vertex
        apex_2_basepointsTetra=Data.apex_2_base(points_Tetra);
        apex_2_base=dot(bar,apex_2_basepointsTetra,2);

        a2bpointsTetra=Data.a2b(points_Tetra);
        a2b=dot(bar,a2bpointsTetra,2);

        r2lpointsTetra=Data.r2l(points_Tetra);
        r2l=dot(bar,r2lpointsTetra,2);

        a2ppointsTetra=Data.a2p(points_Tetra);
        a2p=dot(bar,a2ppointsTetra,2);

        rsinpointsTetra=rsin(points_Tetra);
        rsin=dot(bar,rsinpointsTetra,2);

        rcospointsTetra=rcos(points_Tetra);
        rcos=dot(bar,rcospointsTetra,2);

        r2l_geopointsTetra=Data.r2l_geo(points_Tetra);
        r2l_geo=dot(bar,r2l_geopointsTetra,2);

        Data.a2b_cobi(Data.a2b_cobi==-1)=NaN;
        a2b_cobipointsTetra=Data.a2b_cobi(points_Tetra);
        a2b_cobi=dot(bar,a2b_cobipointsTetra,2);

        a2b_uvcpointsTetra=Data.a2b_uvc(points_Tetra);
        a2b_uvc=dot(bar,a2b_uvcpointsTetra,2);

        lvrv_cobipointsTetra=Data.lvrv_cobi(points_Tetra);
        lvrv_cobi=dot(bar,lvrv_cobipointsTetra,2);

        % a2b_cutpointsTetra=Data.a2b_cut(points_Tetra);
        % a2b_cut=dot(bar,a2b_cutpointsTetra,2);

        F=gradientinterpol(Data.F,length(cube_centroid),points_Tetra,bar);
        F_N=gradientinterpol(Data.F_N,length(cube_centroid),points_Tetra,bar);
        F_S=gradientinterpol(Data.F_S,length(cube_centroid),points_Tetra,bar);

        ahapointsTetra=Data.aha(points_Tetra);
        aha=dot(bar,ahapointsTetra,2);



   %r field
     r=atan2(rsin,rcos)./(2*pi);
     r(rsin<0)=atan2(rsin(rsin<0),rcos(rsin<0))./(2*pi)+1;
     a2b_cobi(isnan(a2b_cobi))=-1;
     r(isnan(r))=-1;

 %write info in points
 
     %save in cells
     cd (monodir)
     save_ensi_MonoAlg_hex_v2(node_vox,elem,Ventricle,d3,Tphi3,tm_cobi,Epiendo,Epiendo3,a2b_uvc,a2b_cobi,r,lvrv_cobi,r2l_geo,a2b,r2l,a2p,F,F_S,F_N,Material,aha);
     
     FastEndo=ones(size(Ventricle));
     FastEndo(Ventricle==1 | Ventricle==-2)=0; %Fast endo is 0. Normal tissue:1
     Healthytissue=Material;%%%%%%%%%%MODIFY IF NECESSARY
     Ik_s=ones(size(Ventricle));  %IKS scaling for personalisation of T wave. Ones by default.

     T=array2table([ALG*10,Tphi3,Epiendo3,a2b,FastEndo,Healthytissue,F,F_S,F_N,Ik_s]);

     writetable(T,strcat(case_name,'.txt'),'WriteVariableNames',0);

     %change format to .alg
     file = fullfile(strcat(case_name,'.txt'));
    [tempDir, tempFile] = fileparts(file); 
    status = copyfile(file, fullfile(tempDir, [tempFile, '.alg']));
    delete (strcat(case_name,'.txt'));
     %create alg file;
    

end

