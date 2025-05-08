% This is script which runs the pipeline from endocardial and epicardial surface meshes to
% simulation files

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
        
Prefix_LV='LV_endo_ex_';   %change according to the format of the input data
Prefix_LV_epi='LV_epi_ex_';   %change according to the format of the input data
Prefix_RV='RV_ex_';   %change according to the format of the input data

meshformat='UKBB'  ; %type of mesh  UKBB--> inputs obtained from the UKBB image data
                   %               cut --> cut biventricular geomtry
                   %               biventricular_open --> biventricular  with open valves
                   %               biventricular_closed  --> biventricular  with closed valves (EM simulations) 

mesh_resolution_fine=1;
mesh_resolution_coarse=1.5;
mesh_resolution_hexa=0.04;
                                 
%% get files name
cd(origpath)
cases = dir([Prefix_LV, '*']);
case_names=string({cases(:).name})';
case_numbers=str2double(extract(case_names,digitsPattern));
total_cases=length(case_numbers);

if ~exist(resultspath,'dir')
  mkdir(resultspath);       
end
cd(resultspath)

%% mesh generation
for index=1:total_cases
        case_number=case_numbers(index);
         
        case_folder = fullfile(resultspath, num2str(case_number));
        input_folder = fullfile(case_folder, 'input');
        %% create case folder of input data
        if ~exist(case_folder, 'dir')
            mkdir(case_folder);
            mkdir(input_folder);
            copyfile(fullfile(origpath, [Prefix_LV, num2str(case_number), '.ply']), input_folder);
            copyfile(fullfile(origpath, [Prefix_LV_epi, num2str(case_number), '.ply']), input_folder);
            copyfile(fullfile(origpath, [Prefix_RV, num2str(case_number), '.ply']), input_folder);
        end
        
        
         %% read mesh

                cd(input_folder);

                %1 create epicardium of right ventricle
                add_RV_width(Prefix_RV,3,case_number); 
                %2 fix LV mitral borders
                fix_LV_borders (Prefix_LV,case_number);
                [original_LV_mesh]=readVTK(fullfile(input_folder, [Prefix_LV, num2str(case_number), '.ply']));%used for labelling

                %3 merge all surfaces
                labels0=merge_surfaces(Prefix_LV,Prefix_LV_epi,Prefix_RV,case_number);
                %4 add lid to LV
                add_lid_LV()
                %5 Remesh surfaces (important step, it fixes mesh issues)
                surf_edge_length=1.5;
                surf0=remesh_surfaces_UKBB(surf_edge_length,case_number);

                %6 Move surfaces meshes to result folder and change directory to
                %that folder
                cd(case_folder)
                movefile(fullfile(input_folder, ['Closed_final_', num2str(case_number), '.ply']), case_folder);
                movefile(fullfile(input_folder, 'labels0.vtk'), case_folder);
           

                %% generate coarse mesh

                %find unit scale
                dist=sqrt(sum((surf0.points-repmat(surf0.points(1,:),length(surf0.points),1)).^2,2));
                scaledist=max(dist); %if ~100 --> mm

                if scaledist >50 &&  scaledist <500
                    disp('scale is in mm ')

                elseif scaledist >5 &&  scaledist <=50
                    disp('scale is in cm -->Units for EM solver and Personalization ')
                    disp('conversion to mm for meshing')
                     surf0.points=surf0.points.*10;
                     labels0.points=labels0.points.*10;

                else
                    error('check your mesh dimensions (not mm nor cm)')
                end
                % if ~exist('coarse.vtu','file')    
                   %coarse tetrahedral  mesh at 1.5 mm resolution.
                   %find inner points of RV and LV
                    %HOLE 1 (LV)
                        labelLV=1;
                        carLV=labels0.cells(labels0.cellData==labelLV,:);
                        pointsLV=labels0.points(unique(carLV),:);
                        pt_LV=mean(pointsLV);                   
                    %HOLE 2 (RV)
                        labelRV=3;
                        carRV=labels0.cells(labels0.cellData==labelRV,:);
                        pointsRV=labels0.points(unique(carRV),:);
                        pt_RV=mean(pointsRV);
                   MeshCoarse=tetrahedral_meshing(surf0,mesh_resolution_coarse,pt_LV,pt_RV);
        
                   disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                   disp('conversion to cm after meshing -->Units for EM solver and Personalization ')
                   MeshCoarse.points=MeshCoarse.points./10;
                   disp('units converted to cm')
                   vtkWrite(MeshCoarse, 'Coarse.vtu');
        
                % else
                %    MeshCoarse=vtkRead('Coarse.vtu');
                % 
                % end

        %% generate fine mesh
        disp('generating fine mesh')
                % if ~exist('fine.vtu','file')
                   MeshFine=tetrahedral_meshing(surf0,mesh_resolution_fine,pt_LV,pt_RV);    
                   disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                   disp('conversion to cm after meshing -->Units for EM solver and Personalization ')
                   MeshFine.points=MeshFine.points./10;
                   surf0.points=surf0.points./10;
                   disp('units converted to cm')
                   vtkWrite(MeshFine, 'fine.vtu');
                % else
                %    MeshFine=vtkRead('fine.vtu');
                % end
       %% check normals 
                disp('Checking normals in volumetric meshes')
                sur_coarse = vtkDataSetSurfaceFilter(MeshCoarse);
                TR_Surf=triangulation(double(sur_coarse.cells),double(sur_coarse.points));
                TR_Vol=triangulation(double(MeshCoarse.cells),double(MeshCoarse.points));
                facenormals = faceNormal(TR_Surf);
                centroid=meshcentroid(double(sur_coarse.points),double(sur_coarse.cells));
               
                % Compute mesh centroid in one face
                epsilon = 1e-4;  % Small offset distance along normal
                offsetPoint = centroid(1,:) + epsilon * facenormals(1,:);
                ID = pointLocation(TR_Vol, offsetPoint);

               
                % Flip if needed
                if ~isnan(ID)
                    disp('Normals are inward-facing, flipping...');
                    MeshCoarse.cells(:,[3 4])=MeshCoarse.cells(:,[4 3]);
                    vtkWrite(MeshCoarse, 'Coarse.vtu');
                    MeshFine.cells(:,[3 4])=MeshFine.cells(:,[4 3]);
                    vtkWrite(MeshFine, 'fine.vtu');
                    sur_coarse = vtkDataSetSurfaceFilter(MeshCoarse);
                else
                    disp('Normals are correctly oriented.');
                end


         %% generate Hexa mesh
          disp('Generating Hexa mesh');
        [ALG,tet_ID,bar]= hexa_mesher('Coarse.vtu', mesh_resolution_hexa,1);
        outputmesh = movefile('hex_Coarse.vtk', ['Hexa_', num2str(case_number), '.vtk']);
        MeshHex = vtkRead(['Hexa_', num2str(case_number), '.vtk']);


        
        %% generate labels
         % if ~exist('labels_final.vtk','file')
            disp('generating labels')
            cd(input_folder)        
            original_LV_mesh.points=original_LV_mesh.points./10;            
            
            opt.original_LV_mesh=original_LV_mesh;
            labelfinal3=Ventricular_Labelling(sur_coarse,meshformat,opt);

            movefile(fullfile(input_folder, 'labels_final.vtk'),fullfile(case_folder, 'labels_final.vtk'));
            % end
               
        
        %% generate fields
        disp('generating fields')
        cd(case_folder)

           %Fiber information
           %Angle Definition
          Fiber_info.Interpolation_in_septum=1; %1--> interpolate setpum points with continuity or disontinuity of fibers
                                                %0--> no interpolation
          Fiber_info.Septum_angleRV=deg2rad(-0);%angle of the septal plane fibers respect to RV
          Fiber_info.Discontinuity_angle=deg2rad(0);%angle between the RV septal fibers and LV septal fibers
          %alpha angle
              Fiber_info.AENDORV=90;
              Fiber_info.AEPIRV=-25;
              Fiber_info.AENDOLV=-60;
              Fiber_info.AEPILV=60;
              Fiber_info.AOTENDOLV=-90;
              Fiber_info.AOTENDORV=90;
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
        Field_generator_UKBB_function24(Fiber_info,meshformat,pericardium_level, epiendo, epiendoRV,case_number);

        cd(case_folder)

        %% creation of generic files for EM simulations
        %this files and formats are created for specific solvers. Creation
        %for a particular solver will required extra formating
        % casepath=strcat(resultspath,'/',num2str(case_number),'\',strcat('ensi',num2str(case_number)));
        % cd(casepath)

        Data4EM=load(fullfile(case_folder, ['ensi_Fine_', num2str(case_number)], 'Case_Fine.mat')) ; %choose between the mesh files  (coarse for EM or fine for EP)

        nameEMFolder=['EM_Control_',num2str(case_number)];
        EMfilespath = fullfile(resultspath, num2str(case_number), nameEMFolder);

        Electrodes_final=electrode_generation(Data4EM,referenceFolder);

        name_EM_case = ['heart_', num2str(case_number)];
        %AlyaFiles_creation_all(EMfilespath,Data4EM,name_EM_case,referenceEM,Electrodes_final,case_number); %functions for Alya not provided
      
        %%generate activation
        mkdir(EMfilespath)
        BCL=0.8;
        generateActivation(name_EM_case,Data4EM,referenceFolder,EMfilespath,'generic',BCL);

        %% Hexa Files generation
%         %read field data
        casepath=fullfile(case_folder, ['ensi', num2str(case_number)]);
        monodir=fullfile(case_folder, 'MonoAlg3D');

        cd(casepath)
        Data=load('Case_coarse.mat');
        %read mesh data

        cd(case_folder)

        mkdir(monodir)
        cd(monodir);
        %generate ALG file
        HexaFieldsGeneration_function_cells_v2(monodir,ALG,MeshCoarse,MeshHex,Data,tet_ID,bar,epiendo,epiendoRV,['Fields_',num2str(case_number)]);


       %new electrodes
       Electrodes_orig=electrode_generation(Data,referenceFolder);
       Electrodes_final_monoAlg=[Electrodes_orig(:,1)-min(Data.v(:,1)),Electrodes_orig(:,2)-min(Data.v(:,2)),Electrodes_orig(:,3)-min(Data.v(:,3))];
       write_vtk_points('Electrodes_final.vtk',Electrodes_final_monoAlg,1:10,1:10);

       %root nodes
       roots_number=4;
       BCL=0.8;
       name_mono=['UKBB_',num2str(case_number)];
       rootnodes=rootnodes_from_IDs(Data,referenceFolder,roots_number);
       rootnodes.time=rootnodes.time-0.02;  %transformation from generic rootnodes (simulations starts at 0.02)
       %conductance
       if mesh_resolution_hexa==0.05
            sigma=[0.000186, 0.000093, 0.000123];          %sigma_l   %sigma_t  %sigma_n

       elseif mesh_resolution_hexa==0.04

             sigma=[0.000310, 0.000155, 0.000205];          %sigma_l   %sigma_t  %sigma_n

       else
           prompt = "input conductance value   ([long, trans, normal]) ";
           sigma = input(prompt);
       end

       %creation of .ini case
       Monoalg_sim=struct;
       Monoalg_sim.BCL=BCL*1000;  %ms
       Monoalg_sim.duration=1000; %ms simulation duration
       Monoalg_sim.num_threads=8;
       Monoalg_sim.dt_pde=0.01;
       Monoalg_sim.output_dir=['./outputs/',name_mono];
       Monoalg_sim.fast_endo_layer_scale=10;
       Monoalg_sim.dt_ode=0.01;
       Monoalg_sim.num_volumes=length(MeshHex.cells);
       Monoalg_sim.original_discretization=mesh_resolution_hexa*10000;
       Monoalg_sim.desired_discretization=mesh_resolution_hexa*10000;

       meshfile=strcat('cluster_mesh_path/',['Fields_',num2str(case_number),'.alg']);

       Monoalg_ini_creation(name_mono,monodir,referenceMonoAlg,meshfile,MeshHex,Monoalg_sim,sigma,Electrodes_final_monoAlg,rootnodes)


       cd(resultspath)
end
