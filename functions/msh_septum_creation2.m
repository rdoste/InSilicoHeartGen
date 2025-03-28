function [face2,label2,label_tetra_msh]=msh_septum_creation2(v2,elem,face,labelp,Tphi,label,body1,label_tetra,body1_Threshold,Surface_tag1,Surface_tag2)
%function that generates the septal surfaces and all the variables required
%for the septum (cobiveco version) field calculation
%body1--> the elements of the RV and LV



              %in order to classify NaN, we look for the neighbors
              TR_mesh=triangulation(double(elem),v2);
              for iter=1:10
                  zeros_body=find(isnan(body1));
                  Neighbors_zeros = neighbors(TR_mesh,zeros_body);
                  Neighbors_val=body1(mode(Neighbors_zeros,2));
                  body1(zeros_body)=Neighbors_val;
              end
              
              
              indxLV=find(body1<body1_Threshold);
              elemLV=elem(indxLV,:);
              surf_septL=volface(elemLV); 
              indxRV=find(body1>body1_Threshold);
              elemRV=elem(indxRV,:);
              surf_septR=volface(elemRV);
              points_sept=intersect(surf_septR,surf_septL);          
              points_RV=find(Tphi==1 | Tphi==-2 ) ;
              points_RV_body=unique(surf_septR);
              points_Sept_RV=setdiff(points_RV_body,points_RV);         
              points_Sept=intersect(points_Sept_RV,points_sept);
              label_sept=zeros(size(v2));
              label_sept(points_Sept)=1;
              septal_surface_label=label_sept(surf_septR);
              septal_surface_index=find(prod(septal_surface_label,2));
              septal_surface=surf_septR(septal_surface_index,:);
              face2= [face;septal_surface];                      
              label2=[label;repmat(Surface_tag1,length(septal_surface),1)]; %          
              label_tetra_msh=ones(size(label_tetra)).*10;
              label_tetra_msh(body1<body1_Threshold)=11;               





end