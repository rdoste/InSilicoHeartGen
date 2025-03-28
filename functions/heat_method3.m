function  [phi2,D,u,X,div_X,phi,pre,B,t]=heat_method3(v2,f2,gamma,DTphi)
%method for reescaling the heat equation results


 %read faces and gradient
 %f2--> faces
 %DTphi-->gradient
 %label_tetra_msh-->labels of tetra
 %Ventricle_label--> 10 is RV, 11 is LV
 
 %gradient in faces
 DTphi_f=zeros(length(f2),3);

 
    for ind=1:length(f2)    
    DTphi_f(ind,:)=DTphi(f2(ind,1),:)+DTphi(f2(ind,2),:)+DTphi(f2(ind,3),:)+DTphi(f2(ind,4),:);
    end
 

 
 
 
%reescale f3 in order to obtain the relative order of the faces
% unique_vertex=unique(f2);
% [~,f4]=ismember(f2,unique_vertex);


 %poisson equation
[D,u,X,div_X,phi,pre,B,t] = heat_geodesic_modified2t(v2,double(f2),gamma,DTphi_f);
   
   %Reescale result between 0 and 1
   
    maxphi = max(phi(:));
    minphi= min(phi(:));
    phi2   = (phi - minphi) / (maxphi - minphi);
 %Result of Poisson equation
 

 
%  
%  %values in tetras
%  phi_elem=phi2(f4);
%  phi_elem=mode(phi_elem,2);

 %write solution
 % write_vtk_rbm('testrv.vtk',N_points,v2,N_faces,f,phi3);  
 % write_vtk_rbm_tetra('testtrv.vtk',length(v3),v3,length(f4),[(1:length(f4))',f4]-1,phi_elem,DTphi_f(label_tetra_msh==Ventricle_label,:)); 
 



end