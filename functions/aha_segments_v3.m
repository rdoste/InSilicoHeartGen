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

function  [AHA]= aha_segments_v3(v2,f2,Ventricle,r,pto,car1,Fid)

        % function that divides the biventricular geometry in the
        % different AHA regions (https://doi.org/10.1161/hc0402.102975)

        % INPUT: 
        % v2: mesh node coordinates
        % f2: tetrahedral mesh connectivity
        % Ventricle: Field with labels from epicardium, miocardium and
        % endocardium of the LV and RV
        % r: circunferential coordinate
        % pto: nodes of the surface mesh
        % car1: faces of the surface mesh
        % Fid: labels of the surface mesh
        
        % OUTPUT:
        % AHA: labels in the nodes of the several RV and LV AHA regions
        


        %Longitudinal direction
        
          v_apex=unique(car1(Fid==5 | Fid==12 ,:));
        %find closest basal point
            pto_base=pto(unique(car1(Fid==4 | Fid==13,:)),:);
            TR_base=delaunayTriangulation(pto_base);
            v_base_id=nearestNeighbor(TR_base,pto(v_apex,:));
            v_base=pto_base(v_base_id,:);            
            center_mitral=mean(v2((Ventricle==-2 ),:));
            center_tricuspid=mean(v2((Ventricle==1 ),:));
            facesbase=car1(Fid==1,:);
            snorm=surfacenorm(pto,facesbase);
            [row,~]=find(isnan(snorm));
            snorm(row,:)=[];
            normal_base=mean(snorm,1);%normal of the basal plane            
            Long_axis_2points=v_base(1,:)-pto(v_apex(1),:);
            if dot( Long_axis_2points,normal_base)<0
             normal_base=normal_base*-1;
            end
               
              
       %%  Left Ventricle 
       %Alignment with z coordinate
            Long_axis=normal_base;
            z=[0 0 1]; 
            b=Long_axis/norm(Long_axis);
            prod_vect=cross(z,b,2);
            angle=atan2(sqrt(sum(prod_vect.^2,2)),dot( z,b,2));%angle between both e1 vectors
            w = cross(z,b);    
            ssc = [0 -w(3) w(2); w(3) 0 -w(1);-w(2) w(1) 0];
            R = eye(3) + ssc + ssc^2*(1-dot(z,b))/(norm(w))^2; 
            
            %alignment of the x axis with the mitral_tricuspid center 
            tri2mit=(center_mitral*R-center_tricuspid*R);
            tri2mit=tri2mit/norm(tri2mit);
            tri2mit=tri2mit(1:2);
            ang=acos(dot(tri2mit,[1 0])); 
            R1=[ cos(ang) -sin(ang) 0; sin(ang) cos(ang) 0; 0 0 1];
            R=R*R1;
            vLV=v2*R;
            zLV=max(vLV(:,3))-min(vLV(:,3));

            LVendo=Ventricle==-2;
                [~,id_max_mitral_endo]=min(vLV(LVendo,3)); 
                p_max_mitral_endo=find(LVendo,id_max_mitral_endo);%endo LV apex point
                p_max_mitral_endo=p_max_mitral_endo(id_max_mitral_endo);
            
            
            RVendo=Ventricle==1;
                [~,id_max_tricusp_endo]=min(vLV(RVendo,3)); 
                p_max_tricusp_endo=find(RVendo,id_max_tricusp_endo);%endo LV apex point
                p_max_tricusp_endo=p_max_tricusp_endo(id_max_tricusp_endo);
                
            RVepi=Ventricle==3;
                [~,id_max_tricusp_epi]=min(vLV(RVepi,3)); 
                p_max_tricusp_epi=find(RVepi,id_max_tricusp_epi);%endo LV apex point
                p_max_tricusp_epi=p_max_tricusp_epi(id_max_tricusp_epi);
                
            
                labelaha=zeros(size(Ventricle));           
                labelaha(vLV(:,3)<(zLV/3+min(vLV(:,3))) & Ventricle<0)=22;
                labelaha(vLV(:,3)>=(zLV/3+min(vLV(:,3))) & Ventricle<0)=33;
                labelaha(vLV(:,3)>=(2*zLV/3+min(vLV(:,3))) & Ventricle<0)=44;
                labelaha((vLV(:,3)-0.02*(abs(max(vLV(:,3))-min(vLV(:,3))))<=vLV(p_max_mitral_endo,3)) & Ventricle<0)=17;
                
                
                labelaha(vLV(:,3)<(zLV/4+min(vLV(:,3))) & Ventricle>0)=55;
                labelaha(vLV(:,3)>=(zLV/4+min(vLV(:,3))) & Ventricle>0)=66;
                labelaha(vLV(:,3)>=(2*zLV/4+min(vLV(:,3))) & Ventricle>0)=77;
                labelaha(vLV(:,3)>=(3*zLV/4+min(vLV(:,3))) & Ventricle>0)=88;


            %% Re-alignment in the LV
            %Alignment with z coordinate
                Long_axis=center_mitral*R-vLV(p_max_mitral_endo,:);
                z=[0 0 1]; 
                b=Long_axis/norm(Long_axis);
                prod_vect=cross(z,b,2);
                angle=atan2(sqrt(sum(prod_vect.^2,2)),dot( z,b,2));%angle between both e1 vectors
                w = cross(z,b);    
                ssc = [0 -w(3) w(2); w(3) 0 -w(1);-w(2) w(1) 0];
                R2 = eye(3) + ssc + ssc^2*(1-dot(z,b))/(norm(w))^2; 
            
            %alignment of the x axis with the mitral_tricuspid center 
                tri2mit=(center_mitral*R*R2-center_tricuspid*R*R2);
                tri2mit=tri2mit/norm(tri2mit);
                tri2mit=tri2mit(1:2);
                ang=acos(dot(tri2mit,[1 0])); 
                R3=[ cos(ang) -sin(ang) 0; sin(ang) cos(ang) 0; 0 0 1];
                R=R*R2*R3;
                vLV=v2*R;
                zLV=max(vLV(:,3))-min(vLV(:,3));


             %% angle division

             %First we calculate the max and min value of angle for the points of
             %the RV
             %After that, we calculate the angle (angle_LV) for the points of the
             %LV in order to obtain the two septal parts

               vect_RV_M=vLV(Ventricle>0,1:2)-repmat(vLV(p_max_mitral_endo,1:2),length(find(Ventricle>0)),1);
               angle_RV_M= atan2(vect_RV_M(:,2),vect_RV_M(:,1));
               angle_RV_M=rad2deg(angle_RV_M);        
               angle_RV_M(angle_RV_M<0)=angle_RV_M(angle_RV_M<0)+360; %angulo de los puntos del RV  

               vect_LV=vLV(:,1:2)-repmat(vLV(p_max_mitral_endo,1:2),length(find(Ventricle)),1);
               angle_LV= atan2(vect_LV(:,2),vect_LV(:,1));
               angle_LV=rad2deg(angle_LV);       
               angle_LV(angle_LV<0)=angle_LV(angle_LV<0)+360;    

               angle_seg=ceil((360-(max(angle_RV_M)-min(angle_RV_M)))/4);%each segment is aproximatly 60�


               %labelling

               labelaha((angle_LV>=min(angle_RV_M) & angle_LV<=(max(angle_RV_M)+min(angle_RV_M))/2) & labelaha==44)=3;
               labelaha((angle_LV<=max(angle_RV_M) & angle_LV>=(max(angle_RV_M)+min(angle_RV_M))/2) & labelaha==44)=2;
               labelaha((angle_LV<=max(angle_RV_M) & angle_LV>=(max(angle_RV_M)+min(angle_RV_M))/2) & labelaha==33)=8;
               labelaha((angle_LV>=min(angle_RV_M) & angle_LV<=(max(angle_RV_M)+min(angle_RV_M))/2) & labelaha==33)=9;
               labelaha((angle_LV<=min(angle_RV_M) & angle_LV>min(angle_RV_M)-angle_seg) & labelaha==44)=4;
               labelaha((angle_LV<=min(angle_RV_M) & angle_LV>min(angle_RV_M)-angle_seg) & labelaha==33)=10;
               labelaha(((angle_LV<=min(angle_RV_M)-angle_seg & angle_LV>=min(angle_RV_M)-angle_seg*2) | angle_LV>=max(angle_RV_M)+angle_seg*2)& labelaha==44)=5;
               labelaha(((angle_LV<=min(angle_RV_M)-angle_seg & angle_LV>=min(angle_RV_M)-angle_seg*2) | angle_LV>=max(angle_RV_M)+angle_seg*2)& labelaha==33)=11;
               labelaha((angle_LV>=max(angle_RV_M) & angle_LV<=max(angle_RV_M)+angle_seg) & labelaha==44)=1;
               labelaha((angle_LV>=max(angle_RV_M) & angle_LV<=max(angle_RV_M)+angle_seg) & labelaha==33)=7;
               labelaha(((angle_LV>=max(angle_RV_M)-angle_seg & angle_LV<=max(angle_RV_M)+angle_seg*2) |  angle_LV<=min(angle_RV_M)-angle_seg*2)& labelaha==44)=6;
               labelaha(((angle_LV>=max(angle_RV_M)-angle_seg & angle_LV<=max(angle_RV_M)+angle_seg*2) |  angle_LV<=min(angle_RV_M)-angle_seg*2) & labelaha==33)=12;
               labelaha((angle_LV<=(max(angle_RV_M)+min(angle_RV_M))/2+45 & angle_LV>=(max(angle_RV_M)+min(angle_RV_M))/2-45) & labelaha==22 )=14;
               labelaha((angle_LV<=(max(angle_RV_M)+min(angle_RV_M))/2+135 & angle_LV>=(max(angle_RV_M)+min(angle_RV_M))/2+45) & labelaha==22 )=13;
               labelaha((angle_LV<=(max(angle_RV_M)+min(angle_RV_M))/2-45 & angle_LV>=(max(angle_RV_M)+min(angle_RV_M))/2-135) & labelaha==22 )=15;
               labelaha((angle_LV>=(max(angle_RV_M)+min(angle_RV_M))/2+135) & labelaha==22 )=16;
               labelaha((angle_LV<=(max(angle_RV_M)+min(angle_RV_M))/2-135) & labelaha==22 )=16;
               
               
               
               
%% Right Ventricle        
 
       %numbers in RV are chosen according with this rule:
       %1 in RV= 17+1=18  ;   2 =19         

       
        if max(Fid)==18
            pto_base_RV=pto(unique(car1(Fid==14,:)),:);
            pto_pulm=pto(unique(car1(Fid==9,:)),:);
            center_Tricuspid=mean(pto_base_RV);
            center_Pulmonary=mean(pto_pulm);
            Long_axis2=center_Pulmonary-center_Tricuspid;
            z2=[0 0 1]; 
            b2=Long_axis2/norm(Long_axis2);
            w2 = cross(z2,b2);    
            ssc2 = [0 -w2(3) w2(2); w2(3) 0 -w2(1);-w2(2) w2(1) 0];
            R2 = eye(3) + ssc2 + ssc2^2*(1-dot(z2,b2))/(norm(w2))^2;
            vRV=v2*R2;
            zRV=max(vRV(:,3))-min(vRV(:,3));           
            Tval=(vRV(:,3)-min(vRV(:,3)))/(max(vRV(:,3))-min(vRV(:,3))); 
    
        end
          
          % labelling 
                   
         labelaha((r>=4/9 & r<5/6) &  labelaha==88)=17+14;
         labelaha(labelaha==88)=17+15;
         if max(Fid)==18
            labelaha(Tval<=0.5 &  (labelaha==17+14 |labelaha==17+15)  )=17+15; %points where circunferential coord is not defined
            labelaha(Tval>0.5 &  (labelaha==17+14 |labelaha==17+15)  )=17+14;  %points where circunferential coord is not defined
         end
         
        labelaha(r>=2/3 & labelaha==77)=17+4;   
        labelaha(r>=2/3 & labelaha==66)=17+9;
        labelaha(r>=2/3 & labelaha==55)=17+13;
        labelaha(r>5/6 & labelaha==17+4)=17+5;
        labelaha(r<=2/9  & labelaha==77)=17+1;
        labelaha((r>2/9 & r<=4/9) & labelaha==77)=17+2;               
        labelaha((r>4/9 & r<=2/3) & labelaha==77)=17+3;
        labelaha((r>2/9 & r<=4/9) & labelaha==66)=17+7;
        labelaha(r<=2/9   & labelaha==66)=17+6;            
        labelaha((r>4/9 & r<=2/3) & labelaha==66)=17+8;                            
        labelaha(r>5/6 & labelaha==17+9)=17+10;                   
        labelaha((r>5/18 & r<=2/3)  &  labelaha==55)=17+12;
        labelaha(r<=5/18 &  labelaha==55)=17+11;


       
       %% Writing     
          
       AHA=labelaha;
           
      % write_vtk_rbm('aha.vtk',length(vLV),vLV,length(f2),f2,AHA);
       

end