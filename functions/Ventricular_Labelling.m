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

function labelfinal3=Ventricular_Labelling(varargin)

        % function that assign labels to different parts of biventricular
        % anatomies

        % INPUT: 
        % sur_coarse: VTK struct of mesh
        % meshtype: type of mesh  
                   %               cut --> biventricular geometry cut in the base by a
                   %               plane
                   %               cut_generic-->  biventricular geometry
                   %               cut in the base (base is not planar)
                   %               open --> biventricular  with open valves
                   %               UKBB  --> closed biventricular geometry with closed valves (from UKBB data) 
        %options:
        % biggestVentRV:         - the biggest ventricle in volume. In healthy conditions is RV (true(1)),  .
                                   % If not, false(1). Most of the "cut" meshes can presente a bigger LV.       
        % original_LV_mesh:      - VTK struct of mesh original surface mesh of the
                                   % LV endocardium, used for mitral valve detection in UKBB meshes

         %anglelid:              -angle used to define lid faces in "cut" geometries
                                   %default (pi/6)
       %  RVseptal_threshold;    - factor used to find endocardial septal RV faces vs wall faces (default 10) 
        
        % OUTPUT:
        % labelfinal3: labels in faces
                    %1 - epicardium
                    %2 - endo LV
                    %3 - endo RV
                    %4 - basal surface (only in cut geometries)
                    %5 - apex (in cut geometries)
                    %6 - septal endocardial wall of the RV
                    %9 - pulmonary valve
                    %10- aortic valve
                    %12- LV apex in epicardium
                    %13- mitral valve
                    %14- tricuspid valve
                    %18- RV apex in endocardium
                    %19 - pulmonary valve in epicardium (in closed meshes)
                    %20- aortic valve in epicardium (in closed meshes)
                    %23- mitral valve in epicardium (in closed meshes)                    
                    %24- tricuspid valve in epicardium (in closed meshes)
                    
                    

        %EXAMPLE OF USAGE:
        %labelfinal3=Ventricular_Labelling(sur_coarse,'UKBB',true(1),original_LV_mesh);


 sur_coarse=varargin{1,1};
 meshtype=string(varargin(2));

% Optional input structure
    if nargin >= 3
        opt = varargin{3};
    else
        opt = struct();
    end

    % Set defaults if fields are missing
    if ~isfield(opt, 'biggestVentRV')
        opt.biggestVentRV = true;
    end

    if ~isfield(opt, 'original_LV_mesh')
        opt.original_LV_mesh = [];
    end

    if ~isfield(opt, 'anglelid')
        opt.anglelid = pi/6;  % angle used to define lid faces in "cut" geometries
    end

    if ~isfield(opt, 'RVseptal_threshold')
        opt.RVseptal_threshold = 10;  %factor used to find endocardial septal RV faces vs wall faces 
    end

    % Check input 
    if nargin > 3
        error('Exceeded the number of function inputs');
    end

    if meshtype == "UKBB" && isempty(opt.original_LV_mesh)
        error('This type of mesh (closed valve mesh) needs the original LV endocardium as input to find the mitral valve');
    end


    %Assign tags to different parts of the ventricle

    node_surf=double(sur_coarse.points);
    face_surf=double(sur_coarse.cells);
    % namemesh='coarse_surface.ply';  %creation of the original surface 
    % save_ply(node_surf,face_surf-1,namemesh(1:end-4));
    %% calculating convexhull 
    disp('CHull calculation')
    Chull = convhull(node_surf);
    write_vtk_surf('Chull.vtk',node_surf,Chull);
    
    %% ray tracing
    
    %edge length
    dist1=sqrt(sum((node_surf(face_surf(:,1),:) - node_surf(face_surf(:,2),:)).^2,2));
    edge_length=median(dist1);
    
    %Remeshing if number of faces is bigger than 100k (speed up
    %calculation)
        if length(face_surf)>100000
        [node_surf_coarse,faces_surf_coarse]=remesh_ventricle(node_surf,face_surf,edge_length*2,1,1,meshtype,1);
        else
          node_surf_coarse=node_surf;
          faces_surf_coarse=face_surf;
        end
   % Ray tracing 
       disp('Ray tracing') 
    
       if meshtype=="cut" || meshtype=="cut_generic"
        vertex1=[node_surf(Chull(:,1),:)];
        vertex2=[node_surf(Chull(:,2),:)];
        vertex3=[node_surf(Chull(:,3),:)];
    
       else
        vertex1=[node_surf_coarse(faces_surf_coarse(:,1),:);node_surf(Chull(:,1),:)];
        vertex2=[node_surf_coarse(faces_surf_coarse(:,2),:);node_surf(Chull(:,2),:)];
        vertex3=[node_surf_coarse(faces_surf_coarse(:,3),:);node_surf(Chull(:,3),:)];
    
       end
    
    
        % tic
          [interp1,interp2,interp3]= calculate_ray_intersections(node_surf,double(face_surf),vertex1,vertex2,vertex3);
        % toc
        % 



        labelf=abs(interp1');
        labelf2=abs(interp2');


          %%

    disp('ventricle labelling...')

    
    %% Assign LV and RV and epicardium surfaces

  if meshtype=="cut" || meshtype=="cut_generic"

    disp('cut mesh ventricle labelling...')
        
        distances=labelf2;
        distances(labelf2==999)=NaN;
        % threshold=max(distances)/10;
        threshold=edge_length;
        
        %threshold calculation                         
        labelfinal=zeros(size(labelf));
        labelfinal(labelf2>threshold & labelf2<999)=1; 
        clusters=numel(finddisconnsurf(face_surf(labelfinal==1,:)));
        cell_lengths=cellfun('length',finddisconnsurf(face_surf(labelfinal==1,:)));
        relsize=max(cell_lengths(cell_lengths<max(cell_lengths)))/max(cell_lengths); %init relsize --> relation between the biggest and the second biggest clusters
        while clusters<2 || relsize<0.3
            threshold=threshold*2;
            labelfinal=zeros(size(labelf));
            labelfinal(labelf2>threshold & labelf2<999)=1;
            clusters=numel(finddisconnsurf(face_surf(labelfinal==1,:)));
            cell_lengths=cellfun('length',finddisconnsurf(face_surf(labelfinal==1,:)));
            relsize=max(cell_lengths(cell_lengths<max(cell_lengths)))/max(cell_lengths);
        end
        

  else
        disp('biventricular mesh ventricle labelling...')
       
        labelfinal=ones(size(interp3'));
        labelfinal(interp3'<=2 & labelf>edge_length*5)=0; % we used labelf to have better results in RV
        
  end
   % write_vtk_surf('test.vtk',node_surf,face_surf,labelfinal); 
    %RV endo vs LV endo
    
    Cluster1=biggestcluster(face_surf(labelfinal==1,:));
    Cluster2=biggestcluster(setdiff(face_surf(labelfinal==1,:),Cluster1,'rows'));
    
    %Bigger Volume
    [~,Volume1] = convhull(node_surf(unique(Cluster1),:));
    [~,Volume2] = convhull(node_surf(unique(Cluster2),:));
    
    
    if Volume2>Volume1
       [ Cluster2,Cluster1] =deal(Cluster1,Cluster2);
    end
    if opt.biggestVentRV==false
        [ Cluster2,Cluster1] =deal(Cluster1,Cluster2);
    end

    %set RV as the one with bigger Convex Hull Volume
    
    labelfinal=ones(size(labelfinal));
    indx_RV0=ismember(face_surf,Cluster1,'rows');
     labelfinal(indx_RV0)=3;
    indx_LV0=ismember(face_surf,Cluster2,'rows');
    labelfinal(indx_LV0)=2;
    
    %% remove isolated patches
    clusters=numel(finddisconnsurf(face_surf(labelfinal==1,:)));
    if clusters>=3
        labelfinal=remove_isolated(labelfinal,node_surf,face_surf,'max');
    end
    
    
    
    %% smoothing
     disp('smoothing labels...')
    %labelfinal=smooth_labelsbi(labelfinal,face_surf,face_surf); %smooth without having the valves
    labelfinal=smooth_labels(labelfinal,node_surf,face_surf); %smooth without having the valves
    labelfinal=remove_isolated(labelfinal,node_surf,face_surf,'max'); 

    % write_vtk_surf('test1.vtk',node_surf,face_surf,labelfinal); 
    
    
    %% Detection of valves
    if meshtype~="UKBB"
        disp('ring detection...')
        labelfinal2 = labelfinal; 
        labelfinal0=labelfinal;
        TR = triangulation(double(face_surf), node_surf);
        centroid = meshcentroid(node_surf, face_surf);
        
        % Identify endocardial and epicardial labels
        LabelEpi = (labelfinal2 == 1);
        LabelEndo = ~LabelEpi; 
        
        % Assign temporary labels
        labelfinal2(LabelEpi) = 10;
        Neighbours = neighbors(TR);
        Label=labelfinal2(Neighbours);
        LabelDiff = abs(labelfinal2(Neighbours) - labelfinal2);
        
        % Detect valve regions
        [ValveRr, ~] = find(LabelDiff == 7);
        [ValveLr, ~] = find(LabelDiff == 8);
        
        % Assign valve neighbors
        NeighborsR = unique(neighbors(TR,ValveRr));
        NeighborsL = unique(neighbors(TR,ValveLr));
        labelfinal2(NeighborsR) = 9;
        labelfinal2(NeighborsL) = 10;
        
        % Restore original epicardial labels
        labelfinal2(LabelEpi) = 1;
        labelfinal3 = labelfinal2;


            %% Define basal plane or valves
              if meshtype=="cut"
            
                
                facenormal = faceNormal(TR);
                facenormal_outer=facenormal(labelfinal==1,:);
                LabelDiff=abs(repmat(labelfinal3,1,3)-Label);
                [ValveRr2,~]=find(LabelDiff==7 );
                [ValveLr2,~]=find(LabelDiff==8);
                NeighborsR2 = unique(neighbors(TR,ValveRr2));
                NeighborsL2 = unique(neighbors(TR,ValveLr2));
                labelfinal3(NeighborsR2)=12;
                labelfinal3(NeighborsL2)=12;
            
                labelfinal3(labelfinal>1 & labelfinal<7)=labelfinal(labelfinal>1 & labelfinal<7);
                %write_vtk_surf('labels_final3.vtk',node_surf,face_surf,labelfinal3);
            
                normal_lid=median(facenormal(labelfinal3==12,:));%we calculate the normal by detecting the valve ring located in the ring and extracting a reference normal
                
                Ref_vector=repmat(normal_lid,length(find(labelfinal==1)),1); 
                prod_vect=cross(Ref_vector,facenormal_outer,2);
                angle=atan2(sqrt(sum(prod_vect.^2,2)),dot( Ref_vector,facenormal_outer,2));%angle between both e1 vectors
                Labelepi_ind=find(LabelEpi);
                
                   %add previously find valves to ensure that endo and epi are never
                    %in contact:
                     labelfinal(NeighborsR)=4;
                     labelfinal(NeighborsL)=4;
                     labelfinal(~LabelEpi)=labelfinal0(~LabelEpi);
            
               
                labelfinal( Labelepi_ind(abs(angle)<opt.anglelid))=4;% angle difference to define lid
            
            
                labelfinal(labelfinal==1)=10;%to avoid issues with the smoothing
                labelfinal=remove_isolated(labelfinal,node_surf,face_surf,'min'); %remove isolated patches
                labelfinal(labelfinal==10)=1;
            
            
              elseif meshtype=="cut_generic"
            
                    %add previously find valves to ensure that endo and epi are never
                    %in contact:
                     labelfinal(NeighborsR)=4;
                     labelfinal(NeighborsL)=4;
                     labelfinal(~LabelEpi)=labelfinal0(~LabelEpi);
            
                 % Define basal plane
                    %creation of a new convex hull with 2 extra points 
                    centerLV=mean(centroid(labelfinal2==2,:));
                    centerLVring=mean(centroid(labelfinal2==10,:));
                    LVaxis=centerLV-centerLVring;
                    ExtrapointLV=centerLVring-LVaxis;
                    
                    centerRV=mean(centroid(labelfinal2==3,:));
                    centerRVring=mean(centroid(labelfinal2==9,:));
                    RVaxis=centerRV-centerRVring;
                    ExtrapointRV=centerRVring-RVaxis;
                    node_surf2=[node_surf;ExtrapointLV;ExtrapointRV]; %new nodes of the second Convex Hull
                    Chull2 = convhull(node_surf2);
                    vertex4=[node_surf2(Chull2(:,1),:)];
                    vertex5=[node_surf2(Chull2(:,2),:)];
                    vertex6=[node_surf2(Chull2(:,3),:)];
                        write_vtk_surf('Chull2.vtk',[node_surf;ExtrapointLV;ExtrapointRV],Chull2);

                    %second raytracing
                    
                    [interp3, interp4]= calculate_ray_intersections(node_surf2,double(face_surf),vertex4,vertex5,vertex6);
                    
                    % labelf3=abs(interp3);
                    % write_vtk_surf('labelsn2.vtk',node_surf,face_surf,labelf2);
                    % write_vtk_surf('labelsn3.vtk',node_surf,face_surf,labelf3);
                    % labelf4=abs(interp4);
                    % write_vtk_surf('labelsn4.vtk',node_surf,face_surf,labelf4);
                    
                    %find differences
                    distancediff=abs(interp2(labelfinal==1)-interp4(labelfinal==1));
                    distancediff2=ones(size(distancediff));
                    distancediff2(distancediff>0)=4;
                    distancediff2(distancediff>500)=1;%points that are outside convexhull 999
                    labelfinal(labelfinal==1)=distancediff2;
             
            
                     
              elseif meshtype=="open"
            
                    %LV valves
                    facecellvalveLV=finddisconnsurf(face_surf(labelfinal2==10,:));
                    [~, max_index] = max(cellfun('size', facecellvalveLV, 1));%detect maximum size array (mitral valve)
                    mitral = cat(1, facecellvalveLV{max_index});
                    facecellvalveLV{:,max_index}=[];%remove the biggest surface
                    [~, max_index] = max(cellfun('size', facecellvalveLV, 1));%detect maximum size array (mitral valve)
                    aortic= cat(1, facecellvalveLV{max_index});
                    %RV valves
                    facecellvalveRV=finddisconnsurf(face_surf(labelfinal2==9,:));
                    [~, max_index] = max(cellfun('size', facecellvalveRV, 1));%detect maximum size array (mitral valve)
                    tricuspid = cat(1, facecellvalveRV{max_index});
                    facecellvalveRV{:,max_index}=[];%remove the biggest surface
                    [~, max_index] = max(cellfun('size', facecellvalveRV, 1));%detect maximum size array (mitral valve)
                    pulmonary= cat(1, facecellvalveRV{max_index});
            
                     [~,col]=ismember(pulmonary,face_surf,'rows');
                     labelfinal2(col)=9;
                     [~,col]=ismember(aortic,face_surf,'rows');
                     labelfinal2(col)=10;
                     [~,col]=ismember(tricuspid,face_surf,'rows');
                     labelfinal2(col)=14;
                     [~,col]=ismember(mitral,face_surf,'rows');
                     labelfinal2(col)=13;
                    
                        
                     pulmonary_v=find(labelfinal2==9);
                     aortic_v=find(labelfinal2==10);
                     tricuspid_v=find(labelfinal2==14);
                     mitral_v=find(labelfinal2==13);
            
            
                       %in order to assure that the valves are classified accordingly
                       % centroid of each valve
                       for inc=[9 10 13 14]
                          center_valve(inc,:)=[(max(centroid(labelfinal2==inc,1))+min(centroid(labelfinal2==inc,1)))/2 ...
                            (max(centroid(labelfinal2==inc,2))+min(centroid(labelfinal2==inc,2)))/2 ...
                            (max(centroid(labelfinal2==inc,3))+min(centroid(labelfinal2==inc,3)))/2];
                       end
                       center_valve( all(~center_valve,2), : ) = [];  
                       Dist = squareform(pdist(center_valve));  
                       %1 and 4--> RV 2 and 3-->LV
                            %obtaining Aortic valve calculating distance
            
                       if Dist(2,1)+Dist(2,4)> Dist(3,1)+Dist(3,4)
                           labelfinal2(aortic_v)=13;
                           labelfinal2(mitral_v)=10;
                           if Dist(2,1)<Dist(2,4)
                               labelfinal2(pulmonary_v)=14;
                               labelfinal2(tricuspid_v)=9;
                           end
                       elseif Dist(3,1)<Dist(3,4)
                               labelfinal2(pulmonary_v)=14;
                               labelfinal2(tricuspid_v)=9;
            
                       end
                       labelfinal=labelfinal2;
            
              end
            
             %% Define RV septal endocardium
                     indx_RV1=(labelfinal==3);
                     indx_RV3=indx_RV1 & labelf<edge_length*opt.RVseptal_threshold; % labelf<edge_length*10;
                     indx_RV6=indx_RV1 & labelf>=edge_length*opt.RVseptal_threshold;% labelf<edge_length*10;
                    
                     Cluster3=biggestcluster(face_surf(indx_RV3,:)); %RV without septum
                     Cluster4=biggestcluster(face_surf(indx_RV6,:)); %RVseptum
                     labelfinal(ismember(face_surf,Cluster3,'rows'))=3;
                     labelfinal(ismember(face_surf,Cluster4,'rows'))=6;

              labelfinals=remove_isolated(labelfinal,node_surf,face_surf,'max');
            
              %% add apex points
             [Fid]=find_Apex(sur_coarse,centroid,labelfinal2);
              if meshtype=="cut" || meshtype=="cut_generic"
                 labelfinals(Fid==12)=5;  % cut geometry
              else
                 labelfinals(Fid==12)=12;  % LV apex
                 labelfinals(Fid==18)=18;  % RV apex
            
              end

            
            
                    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif meshtype=="UKBB"
   %% Redefine RV and LV
     %1-Filter possible epicardial faces
     indx_RV1=(labelfinal==3);
     indx_RV3=indx_RV1 & labelf<edge_length*opt.RVseptal_threshold;
     indx_RV6=indx_RV1 & labelf>edge_length*opt.RVseptal_threshold;

     Cluster3=biggestcluster(face_surf(indx_RV3,:)); %RV without septum
     Cluster4=biggestcluster(face_surf(indx_RV6,:)); %RVseptum
     labelfinal(ismember(face_surf,Cluster3,'rows'))=3;
     labelfinal(ismember(face_surf,Cluster4,'rows'))=6;

    %% Detection of valves
     disp('ring detection...')
     labelfinal2=labelfinal;
     labelfinal2(labelfinal==6)=3;
     TR=triangulation(double(face_surf),node_surf);
     centroid=meshcentroid(node_surf,face_surf);
     labelsurf=sur_coarse;


    %% 
    %LV valves
    %LV mitral valve endo
        labelfinal3=labelfinal;
        node_LV=opt.original_LV_mesh.points;
        face_LV=opt.original_LV_mesh.cells;
        centroid_LV=double(meshcentroid(node_LV,face_LV));
        TR_LV=delaunayTriangulation(centroid_LV);
        [NN,distLV]=nearestNeighbor(TR_LV,centroid);
        Mitral_index0=(labelfinal2==2 & distLV> edge_length);
        Mitral_faces=biggestcluster(face_surf(Mitral_index0,:)); %filter to avoid isolated faces due to LV endo remeshing
        Mitral_index=ismember(face_surf(:,:),Mitral_faces,"rows");

        
        labelfinal3(Mitral_index)=13;
        face_normals=faceNormal(TR);
        %check normals orientation
            % Compute signed volume
            vol = sum(dot(centroid, face_normals, 2)) / 3;
            % Flip if needed
            if vol < 0
                face_normals = -1* face_normals;
            end

     %LV mitral valve epi   (version normals and forcing a circle)
        center_mitral=mean(centroid(labelfinal3==13,:));
        centroid_mitral_endo=meshcentroid(node_surf,face_surf(labelfinal3==13,:));
        distances_mitral=sqrt(sum((centroid_mitral_endo-repmat(center_mitral,length(centroid_mitral_endo),1)).^2,2));
        radius_mitral= prctile(distances_mitral,95); % we don't take the maximum value to avoid any outliers to interfere
        %project values to epi
        centroid_epi=meshcentroid(node_surf,face_surf(labelfinal3==1,:));
        TR_epi=delaunayTriangulation(centroid_epi);
        [center_mitral_epi_index]=nearestNeighbor(TR_epi,center_mitral);
        center_mitral_epi=centroid_epi(center_mitral_epi_index,:);
        distances_mitral_epi=sqrt(sum((centroid_epi-repmat(center_mitral_epi,length(centroid_epi),1)).^2,2));
        label_dist_epi=ones(1,length(distances_mitral_epi));
        label_dist_epi(distances_mitral_epi<radius_mitral)=23;
        labelfinal3(labelfinal3==1)=label_dist_epi;


    %LV apex
warning('off','all')
        LVaxis=median(face_normals(labelfinal3==13,:));
        surfLV.cells=face_surf(labelfinal2==2,:);
        surfLV.points=node_surf;
        surfLV.cellTypes=uint8(ones(length(face_surf(labelfinal2==2)),1)*5);
        longAxLV = computeLongAxis(surfLV, LVaxis);
   %project values along long axis

        labelsurf.pointData.heightLV = double(labelsurf.points*longAxLV');
        labelsurf= vtkPointDataToCellData(labelsurf);
        % vtkWrite(labelsurf, 'labelsurf.vtp');

     %Find maximum value in lV and RV to define Apex

        [~,ApexLV_id_LV]=max(labelsurf.cellData.heightLV(labelfinal2==1));
        ApexLV_id=find(labelfinal2==1);
        ApexLV_id=ApexLV_id(ApexLV_id_LV);
        labelfinal3(ApexLV_id)=12;




    %RV apex
    
        RVaxis=median(face_normals(labelfinal3==13,:));
        surfRV.cells=face_surf(labelfinal2==3,:);
        surfRV.points=node_surf;
        surfRV.cellTypes=uint8(ones(length(face_surf(labelfinal2==3)),1)*5);
        longAxRV = computeLongAxis(surfRV, RVaxis);
   %project values along long axis

        labelsurf.pointData.heightRV = double(labelsurf.points*longAxRV');
        labelsurf= vtkPointDataToCellData(labelsurf);
        % vtkWrite(labelsurf, 'labelsurf2.vtp');

     %Find maximum value in lV and RV to define Apex

       [~,ApexRV_id_RV]=max(labelsurf.cellData.heightRV(labelfinal2==3));
        ApexRV_id=find(labelfinal2==3);
        ApexRV_id=ApexRV_id(ApexRV_id_RV);
        labelfinal3(ApexRV_id)=18;
 
%%  find pulmonary 
    LV2RV=mean(centroid(labelfinal3==3,:))-mean(centroid(labelfinal3==2,:));
    
     %find anterior 2 posterior direction
      RVaxispulm=1*cross(LV2RV,RVaxis);
    % RVaxispulm=LV2RV;
        surfRV.cells=face_surf(labelfinal2==3,:);
        surfRV.points=node_surf;
        surfRV.cellTypes=uint8(ones(length(face_surf(labelfinal2==3)),1)*5);
        RVaxispulm = computeLongAxis(surfRV, RVaxispulm);
       %project values along long axis
        labelsurf3=labelsurf;
        labelsurf3.pointData.heightRV2 = double(labelsurf.points*RVaxispulm');
         labelsurf3= vtkPointDataToCellData(labelsurf3);
             % vtkWrite(labelsurf3, 'labelsurf3.vtp');
    warning('on','all')
        %find center of pulmonary and calculate geodesic
        [~,Pulm_id_RV22]=max(labelsurf3.cellData.heightRV2(labelfinal2==3));
        PulmRV_id=find(labelfinal2==3);
        PulmRV_id=PulmRV_id(Pulm_id_RV22);
          Pulm_RV_point=face_surf(PulmRV_id,1);
      [~,Pulm_id_sept]=max(labelsurf3.cellData.heightRV2(labelfinal3==6));
        PulmRV_id_sept=find(labelfinal3==6);
        PulmRV_id_sept=PulmRV_id_sept(Pulm_id_sept);
        Pulm_RV_sept_point=face_surf(PulmRV_id_sept,1);
    
        %find closest point, between both points
       [Pulm_center_id,~]=nearestNeighbor(TR,(node_surf(Pulm_RV_point,:)+node_surf(Pulm_RV_sept_point,:))./2);
    
      [~,~,~,~,phipulm,~] = heat_geodesic(node_surf,face_surf,Pulm_center_id,[],'BoundaryConditions','robin','Legacy',1);
      phipulm_norm=(phipulm-min(phipulm))./(max(phipulm)-min(phipulm)); 
    
      labelsurf.pointData.phipulm_norm=phipulm_norm;
      labelsurf= vtkPointDataToCellData(labelsurf);
      % define maximum value of pulmonary valve
      max_limit_pulm=min(labelsurf.cellData.phipulm_norm(labelfinal3==6));
      labelfinal3(labelsurf.cellData.phipulm_norm<max_limit_pulm & labelfinal3==3 )=9;
    
    
%%
% pulmary valve epi
    centroid_pulm_endo=meshcentroid(node_surf,face_surf(labelfinal3==9,:));
        TR_pulm_endo=delaunayTriangulation(centroid_pulm_endo);
        [~,dist_pulm]=nearestNeighbor(TR_pulm_endo,centroid);
        labelfinal3(labelfinal2==1 & dist_pulm< 0.35)=19;
    
    
 %% detect tricuspid valve endo, using the plane of the mitral surface
      labelsurf2=labelsurf;
      labelsurf2.pointData.heightLV = double(labelsurf2.points*LVaxis');
      labelsurf2= vtkPointDataToCellData(labelsurf2);
               % vtkWrite(labelsurf2, 'labelsurf4.vtp');

%%
       for percentages=10:-0.5:1
           label_aux=labelfinal3;
           P2 = prctile(labelsurf2.cellData.heightLV(labelfinal2==3),percentages);
           label_aux(labelfinal2==3 & labelsurf2.cellData.heightLV <= P2 )=14;
           number_clusters_Tricuspid=numel(finddisconnsurf(face_surf(label_aux==14,:)));
           if number_clusters_Tricuspid>1
               break
           end
       end
        clusters_Tricuspid=(finddisconnsurf(face_surf(label_aux==14,:)));
    
           cluster_Tricusp_1= cat(1, clusters_Tricuspid{1});
           [~,col_index1]=ismember(cluster_Tricusp_1,face_surf,'rows');
           values_cluster1=mean(labelsurf.cellData.phipulm_norm(col_index1)); %we discern between two clusters using the distance from the LV apex
           index_Tricusp=col_index1;
    
         if number_clusters_Tricuspid>1 
           cluster_Tricusp_2= cat(1, clusters_Tricuspid{2});
           [~,col_index2]=ismember(cluster_Tricusp_2,face_surf,'rows');
           values_cluster2=mean(labelsurf.cellData.phipulm_norm(col_index2)); %we discern between two clusters using the distance from the LV apex 
           if values_cluster2 > values_cluster1
                    index_Tricusp=col_index2;
           end
         end
    
        labelfinal3(index_Tricusp)=14;


    % tricuspid valve epi
        centroid_tri_endo=meshcentroid(node_surf,face_surf(labelfinal3==14,:));
        TR_tri_endo=delaunayTriangulation(centroid_tri_endo);
        [~,dist_tri]=nearestNeighbor(TR_tri_endo,centroid);
        labelfinal3(labelfinal2==1 & dist_tri< 0.35)=24;

        labelfinals=labelfinal3;

    end

  %%

    write_vtk_surf('labels_final.vtk',node_surf,face_surf,labelfinals);


    %functions
        function        C=biggestcluster(MatrixConnectivity)
                    %use connectivity matrix for getting the biggest cluster of faces
                      clusters=finddisconnsurf(MatrixConnectivity);              
                      [~, max_index] = max(cellfun('size', clusters, 1));%detect maximum size array 
                      C= cat(1, clusters{max_index});
    
        end
        function [Fid2,labelsurf]=find_Apex(labelsurf,centroid,Fid)
                 %find apices from left and right ventricle 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Compute Long Axis 
                Fid2=Fid;
            
                centerLV=mean(centroid(Fid==2,:));
                centerLVring=mean(centroid(Fid==10 | Fid==13,:));
                LVaxis=centerLV-centerLVring;
                
                centerRV=mean(centroid(Fid==3,:));
                centerRVring=mean(centroid(Fid==9 | Fid==14,:));
                RVaxis=centerRV-centerRVring;
                
                surfLV.cells=labelsurf.cells(Fid==2,:);
                surfLV.points=labelsurf.points;
                surfLV.cellTypes=uint8(ones(length(labelsurf.cells(Fid==2)),1)*5);
                surfRV.cells=labelsurf.cells(Fid==3 | Fid==6,:);
                surfRV.points=labelsurf.points;
                surfRV.cellTypes=uint8(ones(length(labelsurf.cells(Fid==3 | Fid==6)),1)*5);
                
                longAxLV = computeLongAxis(surfLV, LVaxis);
                longAxRV = computeLongAxis(surfRV, RVaxis);
            
               %project values along long axis
            
                labelsurf.pointData.heightLV = double(labelsurf.points*longAxLV');
                labelsurf.pointData.heightRV = double(labelsurf.points*longAxRV');
                
                labelsurf= vtkPointDataToCellData(labelsurf);
            
            %Find maximum value in lV and RV to define Apex
            
                [~,ApexLV_id_LV]=max(labelsurf.cellData.heightLV(Fid==1));
                ApexLV_id=find(Fid==1);
                ApexLV_id=ApexLV_id(ApexLV_id_LV);
                
                [~,ApexRV_id_RV]=max(labelsurf.cellData.heightRV(Fid==3 | Fid==6));
                ApexRV_id=find(Fid==3 | Fid==6);
                ApexRV_id=ApexRV_id(ApexRV_id_RV);

            %Assign apex values
            
                Fid2(ApexLV_id)=12;
                Fid2(ApexRV_id)=18;
        
    
    end

end

