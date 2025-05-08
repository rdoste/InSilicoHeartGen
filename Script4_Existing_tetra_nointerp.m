% This is script which runs the pipeline from an existing cut biventricular tetrahedral mesh 

clear
%% add folders with matlab libraries and functions
%set one drive path
current_path=pwd;
addpath(genpath(fullfile(current_path, 'functions')));
addpath(genpath(fullfile(current_path, 'dependencies')));

%% configure input data

origpath = fullfile(current_path, 'inputs'); % path with the surface meshes
resultspath = fullfile(current_path, 'outputs'); % path with the resulting files and fields
referenceFolder = fullfile(current_path, 'functions', 'reference_activation'); %reference folder for activation and electrodes
referenceMonoAlg= fullfile(current_path, 'functions', 'reference_activation'); %reference folder for ini files generation
        

meshformat='cut'  ; %type of mesh  UKBB--> inputs obtained from the UKBB image data
                   %                cut --> cut biventricular geomtry
                   %                open --> biventricular  with open valves
                   %                closed  --> biventricular  with closed valves (EM simulations) 

                                 
%% get files name
cd(origpath)
name_origin='cut_mesh_tetra.vtu'; %name of the original surface
name_final='Patient_3';%name of the final mesh and folder
if ~exist(resultspath,'dir')
  mkdir(resultspath);       
end
cd(resultspath)

%% mesh generation
for index=1

    %% read original mesh

                MeshCoarse=vtkRead(fullfile(origpath,name_origin));        
                sur_coarse = vtkDataSetSurfaceFilter(MeshCoarse);

                case_folder = fullfile(resultspath, name_final);
                mkdir(case_folder)
                cd(case_folder)


                %% generate coarse mesh

                %find unit scale
                dist=sqrt(sum((sur_coarse.points-repmat(sur_coarse.points(1,:),length(sur_coarse.points),1)).^2,2));
                scaledist=max(dist); %if ~100 --> mm

                if scaledist >50 &&  scaledist <500
                    disp('scale is in mm ')

                elseif scaledist >5 &&  scaledist <=50
                    disp('scale is in cm -->Units for EM solver and Personalization ')
                    disp('conversion to mm for meshing')
                     sur_coarse.points=sur_coarse.points.*10;
                     MeshCoarse.points=MeshCoarse.points.*10;

                else
                    error('check your mesh dimensions (not mm nor cm)')
                end
              
      
        %% check normals 
                disp('Checking normals in volumetric meshes')
                TR_Surf=triangulation(double(sur_coarse.cells),double(sur_coarse.points));
                TR_Vol=triangulation(double(MeshCoarse.cells),double(MeshCoarse.points));
                facenormals = faceNormal(TR_Surf);
                centroid=meshcentroid(double(sur_coarse.points),double(sur_coarse.cells));
               
                % Compute mesh centroid in one face
                epsilon = 1e-4;  % Small offset distance along normal
                offsetPoint = facenormals(1,:) + epsilon * facenormals(1,:);
                ID = pointLocation(TR_Vol, offsetPoint);

               
                % Flip if needed
                if isnan(ID)
                    disp('Normals are inward-facing, flipping...');
                    MeshCoarse.cells(:,[3 4])=MeshCoarse.cells(:,[4 3]);
                    vtkWrite(MeshCoarse, 'Coarse.vtu');
                    sur_coarse = vtkDataSetSurfaceFilter(MeshCoarse);
                else
                    disp('Normals are correctly oriented.');
                    vtkWrite(MeshCoarse, 'Coarse.vtu');
                end



      
        
        %% generate labels
         % if ~exist('labels_final.vtk','file')
            disp('generating labels')
            opt.biggestVentRV=true(1); %the biggest ventricle is RV (true(1). If not, false(1). Most of the "cut" meshes can presente a bigger LV. 
            opt.RVseptal_threshold=7; %factor used to find endocardial septal RV faces vs wall faces 
            labelfinal3=Ventricular_Labelling(sur_coarse,meshformat,opt);

         % end
            %TO DO: add check of neighbors of label 6 to force repeat with
            %biggestVentRV changed

        
        
        %% generate fields
        cd(case_folder)
        disp('generating fields')
           %Fiber information
           %Angle Definition
          Fiber_info.Interpolation_in_septum=1; %1--> interpolate setpum points with continuity or disontinuity of fibers
                                                %0--> no interpolation
          Fiber_info.Septum_angleRV=deg2rad(-0);%angle of the septal plane fibers respect to RV
          Fiber_info.Discontinuity_angle=deg2rad(0);%angle between the RV septal fibers and LV septal fibers
          %alpha angle
              Fiber_info.AENDORV=80;
              Fiber_info.AEPIRV=-40;
              Fiber_info.AENDOLV=-60;
              Fiber_info.AEPILV=60;
              Fiber_info.AOTENDOLV=0;
              Fiber_info.AOTENDORV=0;
              Fiber_info.AOTEPILV=0;
              Fiber_info.AOTEPIRV=0;
          %beta angle
           Fiber_info.beta=1; % 1--> with transverse angle    0--> without angle    
                  Fiber_info.beta_endoRV =deg2rad(180); 
                  Fiber_info.beta_epiRV =deg2rad(180);  
                  Fiber_info.beta_endoLV =deg2rad(180);
                  Fiber_info.beta_epiLV=20;   %difference between epi/endo and the center of the LV


        %generate fields
        pericardium_level=0.8;
        epiendo=[70 0 30]; % percentage of endo/ mid/ epi#
        epiendoRV=[70 0 30]; % percentage of endo/ mid/ epi (RV septal wall as epi)
        requiresInterpolation=0; %0--> no need interpolation for a large number of points
                                 %1--> interpolation requiered    
        Field_generator_UKBB_function24(Fiber_info,meshformat,pericardium_level, epiendo, epiendoRV,requiresInterpolation,[]);

        cd(case_folder)

        
       cd(resultspath)
end

