function [a2b,r2l,a2p,a2b_vector,r2l_vector,a2p_vector]=ventricular_directions(v,pto,car,Fid,face2,labelf2)
    %function that calculate 3 different gradients in the ventricle
    %1-Apexbase
    % 2-RV to LV
    % 3-Anterior to posterior
    car1=car+1;
    %a2b
    %find basal plane
    if max(Fid)==5 || max(Fid)==6
            facesbase=car1(Fid==4,:);
            snorm=surfacenorm(pto,facesbase);
            normal_base=mean(snorm,1);%normal of the basal plane
            %Rotation to obtain the Apex to base values projected in the z axis
            Long_axis=normal_base;
            z2=[0 0 1]; 
            b2=Long_axis/norm(Long_axis);
            w2 = cross(z2,b2);    
            ssc2 = [0 -w2(3) w2(2); w2(3) 0 -w2(1);-w2(2) w2(1) 0];
            R2 = eye(3) + ssc2 + ssc2^2*(1-dot(z2,b2))/(norm(w2))^2;
            vRV=v*R2;         
            a2b=(vRV(:,3)-min(vRV(:,3)))/(max(vRV(:,3))-min(vRV(:,3))); 
    elseif max(Fid)==18
          pto_base=pto(unique(car1(Fid==13,:)),:);
          pto_apex=pto(unique(car1(Fid==12,:)),:);
          center_apex=mean(pto_apex);
          center_base=mean(pto_base);
          Long_axis=center_base-center_apex;
          z2=[0 0 1]; 
          b2=Long_axis/norm(Long_axis);
          w2 = cross(z2,b2);    
          ssc2 = [0 -w2(3) w2(2); w2(3) 0 -w2(1);-w2(2) w2(1) 0];
          R2 = eye(3) + ssc2 + ssc2^2*(1-dot(z2,b2))/(norm(w2))^2;
          vRV=v*R2;         
          a2b=(vRV(:,3)-min(vRV(:,3)))/(max(vRV(:,3))-min(vRV(:,3))); 
    elseif max(Fid)==24 %a2b in the direction of the Lid normal
          TR=triangulation(double(car1),pto);
          face_normals=faceNormal(TR);
          Long_axis=-median(face_normals(Fid==13,:));
          z2=[0 0 1]; 
          b2=Long_axis/norm(Long_axis);
          w2 = cross(z2,b2);    
          ssc2 = [0 -w2(3) w2(2); w2(3) 0 -w2(1);-w2(2) w2(1) 0];
          R2 = eye(3) + ssc2 + ssc2^2*(1-dot(z2,b2))/(norm(w2))^2;
          vRV=v*R2;         
          a2b=(vRV(:,3)-min(vRV(:,3)))/(max(vRV(:,3))-min(vRV(:,3))); 
    
    end
    a2b_vector=b2;
    
    %r2l
    facesbase2=face2(labelf2==20,:);
    snorm2=surfacenorm(v,facesbase2);
    normal_septum=mean(snorm2,1);%normal of septum
    Long_axis2=normal_septum;
    z2=[0 0 1]; 
    b2=Long_axis2/norm(Long_axis2);
    w2 = cross(z2,b2);    
    ssc2 = [0 -w2(3) w2(2); w2(3) 0 -w2(1);-w2(2) w2(1) 0];
    R2 = eye(3) + ssc2 + ssc2^2*(1-dot(z2,b2))/(norm(w2))^2;
    vRV=v*R2;         
    r2l=(vRV(:,3)-min(vRV(:,3)))/(max(vRV(:,3))-min(vRV(:,3))); 
    
    
    r2l_vector=b2;
    
    
    %a2p
    
    anterior=cross(Long_axis,normal_septum);
    
    Long_axis2=anterior;
    z2=[0 0 1]; 
    b2=Long_axis2/norm(Long_axis2);
    w2 = cross(z2,b2);    
    ssc2 = [0 -w2(3) w2(2); w2(3) 0 -w2(1);-w2(2) w2(1) 0];
    R2 = eye(3) + ssc2 + ssc2^2*(1-dot(z2,b2))/(norm(w2))^2;
    vRV=v*R2;         
    a2p=(vRV(:,3)-min(vRV(:,3)))/(max(vRV(:,3))-min(vRV(:,3))); 
    a2p_vector=b2;
end