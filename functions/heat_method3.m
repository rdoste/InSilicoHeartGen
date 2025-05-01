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

 %poisson equation
[D,u,X,div_X,phi,pre,B,t] = heat_geodesic_modified2t(v2,double(f2),gamma,DTphi_f);
   
   %Reescale result between 0 and 1
   
    maxphi = max(phi(:));
    minphi= min(phi(:));
    phi2   = (phi - minphi) / (maxphi - minphi);
 

%  %values in tetras
%  phi_elem=phi2(f4);
%  phi_elem=mode(phi_elem,2);



end