function labelf=smooth_labels(label,node_surf,face_surf)
        labelf=label;  

        TR=triangulation(double(face_surf),node_surf);

        NE=neighbors(TR); %neighbors of each triangle
        %in case of NaNs in NE:
        [row, col] = find(isnan(NE));
        in=sub2ind(size(NE), row,col);
        NE(in)=NE(in+1);

             
%% smooth labels using nearest neighbors values
        NE_2=[  neighbors(TR,NE(:,1)) neighbors(TR,NE(:,2)) neighbors(TR,NE(:,3))];
        
        Lab_NE=labelf(NE);
        Lab_NE_2=labelf(NE_2);
        SMF=[Lab_NE,Lab_NE_2];
        [M2,F2]= mode(SMF,2);
        % indF=find(F2<8);
        % labelf(indF)=M2(indF);


        ind=0;

         while ind <50
           Lab_NE=M2(NE);
           Lab_NE_2=M2(NE_2);
           SMF=[Lab_NE,Lab_NE_2]; 
           [M2,F2]= mode(SMF,2);
           ind=ind+1;
         end

         labelf=M2;



       

  % 
   % % remove disconected patches
   %      label_number=unique(label);
   % 
   % 
   % 
   % 
   %                      for i=1:numel(label_number)
   %                          label_ind=label_number(i);
   % 
   %                          facecell=finddisconnsurf(face_surf(labelf==label_ind,:));%find all faces 
   %                          while size(facecell,2)>1
   %                              [~, max_index] = max(cellfun('size', facecell, 1));%detect maximum size array
   %                              facecell{:,max_index}=[];%remove the biggest surface
   % 
   %                              labelsf_isolated =ismember(face_surf,cat(1, facecell{:}),'rows');%isolated faces
   %                              NE_2=[  neighbors(TR,find(labelsf_isolated))];
   %                              NE2_2=[  neighbors(TR,NE_2(:,1)) neighbors(TR,NE_2(:,2)) neighbors(TR,NE_2(:,3))];
   % 
   % 
   %                              SM3=labelf(NE2_2);
   %                              SM4=SM3;
   %                              SM4(SM3==label_ind)=NaN; %remove current label
   %                              M3= (max(SM4'))'; %calculate the mode of the neighbors and substitute
   %                             labelf(labelsf_isolated)=M3;
   %                             facecell=finddisconnsurf(face_surf(labelf==label_ind,:));%find all faces 
   % 
   %                          end
   %                      end

       
end