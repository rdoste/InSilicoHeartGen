function labelf=remove_isolated(label,node_surf,face_surf,method)
%method max --> substitutes isolated patches using the neighbors max values
%method min --> substitutes isolated patches using the neighbors min values
%method mod --> substitutes isolated patches using the neighbors mode values

   labelf=label;  
   TR=triangulation(double(face_surf),node_surf);

   % remove disconected patches
        label_number=unique(label);

switch method
    case 'max'
                        for i=1:numel(label_number)
                            label_ind=label_number(i);
                            labelf(labelf==label_ind)=0; %force the current label to be the minimum
                            facecell=finddisconnsurf(face_surf(labelf==0,:));%find all faces 
                            iter=0;
                            while size(facecell,2)>1 && iter<20
                                [~, max_index] = max(cellfun('size', facecell, 1));%detect maximum size array
                                facecell{:,max_index}=[];%remove the biggest surface

                                labelsf_isolated =ismember(face_surf,cat(1, facecell{:}),'rows');%isolated faces
                                NE2=[  neighbors(TR,find(labelsf_isolated))];

                                SM2=labelf(NE2);
                                M2= (max(SM2'))'; %calculate the max of the neighbors and substitute
                               labelf(labelsf_isolated)=max(M2);
                               facecell=finddisconnsurf(face_surf(labelf==0,:));%find all faces 
                               iter=iter+1;

                            end
                            labelf(labelf==0)=label_ind;
                        end

    case 'min'
                        for i=1:numel(label_number)
                            label_ind=label_number(i);

                            facecell=finddisconnsurf(face_surf(labelf==label_ind,:));%find all faces 
                            iter=0;
                            while size(facecell,2)>1 && iter<20
                                [~, max_index] = max(cellfun('size', facecell, 1));%detect maximum size array
                                facecell{:,max_index}=[];%remove the biggest surface

                                labelsf_isolated =ismember(face_surf,cat(1, facecell{:}),'rows');%isolated faces
                                NE2=[  neighbors(TR,find(labelsf_isolated))];

                                SM2=labelf(NE2);
                                M2= ( min(SM2'))'; %calculate the min of the neighbors and substitute
                               labelf(labelsf_isolated)=min(M2);
                               facecell=finddisconnsurf(face_surf(labelf==label_ind,:));%find all faces 
                               iter=iter+1;

                            end
                        end


    case 'mod'

                        for i=1:numel(label_number)
                            label_ind=label_number(i);

                            facecell=finddisconnsurf(face_surf(labelf==label_ind,:));%find all faces 
                            iter=0;
                            while size(facecell,2)>1 && iter<20
                                [~, max_index] = max(cellfun('size', facecell, 1));%detect maximum size array
                                facecell{:,max_index}=[];%remove the biggest surface

                                labelsf_isolated =ismember(face_surf,cat(1, facecell{:}),'rows');%isolated faces
                                NE2=[  neighbors(TR,find(labelsf_isolated))];

                                SM2=labelf(NE2);
                                M2= ( mode(SM2'))'; %calculate the min of the neighbors and substitute
                               labelf(labelsf_isolated)=M2;
                               facecell=finddisconnsurf(face_surf(labelf==label_ind,:));%find all faces 
                               iter=iter+1;

                            end
                        end

    otherwise
        error('Unexpected method')
       
       
end