% This is a sample script which runs the pipeline from an existing cut biventricular surface mesh to
% simulation files

clear
%% add folders with matlab libraries and functions
%set one drive path
current_path=pwd;
addpath (genpath(strcat(current_path,'\functions')));
addpath (genpath(strcat(current_path,'\dependencies')));

%% configure input data

origpath=strcat(current_path,'\inputs\'); % path with the surface meshes
resultspath=strcat(current_path,'\outputs\'); % path with the resulting files and fields
referenceFolder=strcat(current_path,'\functions\reference_activation'); %reference folder for activation and electrodes
referenceMonoAlg=strcat(current_path,'\functions\reference_activation'); %reference folder for ini files generation
        

meshformat='cut'  ; %type of mesh  UKBB--> inputs obtained from the UKBB image data
                   %                cut --> cut biventricular geomtry
                   %                open --> biventricular  with open valves
                   %                closed  --> biventricular  with closed valves (EM simulations) 

mesh_resolution_fine=1;
mesh_resolution_coarse=1.5;
mesh_resolution_hexa=0.04;
                                 
%% get files name
cd(origpath)
name_origin='cut_mesh.ply'; %name of the original surface
name_final='Patient_1';%name of the final mesh and folder
if ~exist(resultspath,'dir')
  mkdir(resultspath);       
end
cd(resultspath)

%% mesh generation
for index=1

    %% read original mesh

                surf0=vtkRead(strcat(origpath,name_origin));        

                mkdir(strcat(resultspath,name_final))
                cd(strcat(resultspath,name_final))


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

                else
                    error('check your mesh dimensions (not mm nor cm)')
                end
                % if ~exist('coarse.vtu','file')    
                   %coarse tetrahedral 
                   MeshCoarse=tetrahedral_meshing(surf0,mesh_resolution_coarse,[],[]);
        
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
                   MeshFine=tetrahedral_meshing(surf0,mesh_resolution_fine,[],[]);    
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
                offsetPoint = facenormals(1,:) + epsilon * facenormals(1,:);
                ID = pointLocation(TR_Vol, offsetPoint);

               
                % Flip if needed
                if isnan(ID)
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
          
        [ALG,tet_ID,bar]= hexa_mesher('Coarse.vtu', mesh_resolution_hexa,1);
        outputmesh=movefile('hex_Coarse.vtk',strcat('Hexa_',name_final,'.vtk'));
        MeshHex=vtkRead(strcat('Hexa_',name_final,'.vtk'));



        
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
        cd(strcat(resultspath,name_final))

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
        Field_generator_UKBB_function24(Fiber_info,meshformat,pericardium_level, epiendo, epiendoRV);

        cd(strcat(resultspath,name_final))

        %% creation of generic files for EM simulations
        %this files and formats are created for specific solvers. Creation
        %for a particular solver will required extra formating
        % casepath=strcat(resultspath,'/',name_final,'\',strcat('ensi',name_final));
        % cd(casepath)
        Data4EM=load(strcat(resultspath,name_final,'\ensi_Fine_','\Case_Fine.mat')) ; %choose between the mesh files  (coarse for EM or fine for EP)

        nameEMFolder=strcat('EM_Control_',name_final);
        EMfilespath=strcat(resultspath,name_final,'\',nameEMFolder);

        Electrodes_final=electrode_generation(Data4EM,referenceFolder);

        name_EM_case=strcat('heart_',name_final);
        %AlyaFiles_creation_all(EMfilespath,Data4EM,name_EM_case,referenceEM,Electrodes_final,case_number); %functions for Alya not provided

        %%generate activation
        mkdir(EMfilespath)
        BCL=0.8;
        generateActivation(name_EM_case,Data4EM,referenceFolder,EMfilespath,'generic',BCL);



        %% Hexa Files generation
%         %read field data
        casepath=strcat(resultspath,name_final,'\',strcat('ensi','\'));
        monodir=strcat(resultspath,name_final,'\MonoAlg3D');

        cd(casepath)
        Data=load('Case_coarse.mat');
        %read mesh data

       cd(strcat(resultspath,name_final))

        mkdir(monodir)
        cd(monodir);
        %generate ALG file
        HexaFieldsGeneration_function_cells_v2(monodir,ALG,MeshCoarse,MeshHex,Data,tet_ID,bar,epiendo,epiendoRV);


        %new electrodes
        Electrodes_orig=electrode_generation(Data,referenceFolder);
        Electrodes_final_monoAlg=[Electrodes_orig(:,1)-min(Data.v(:,1)),Electrodes_orig(:,2)-min(Data.v(:,2)),Electrodes_orig(:,3)-min(Data.v(:,3))];
        write_vtk_points('Electrodes_final.vtk',Electrodes_final_monoAlg,1:10,1:10);

       %root nodes
       roots_number=4;
       BCL=0.8;
       name_mono=strcat('UKBB_',name_final);
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
       Monoalg_sim.output_dir=strcat('./outputs/',name_mono);
       Monoalg_sim.fast_endo_layer_scale=10;
       Monoalg_sim.dt_ode=0.01;
       Monoalg_sim.num_volumes=length(MeshHex.cells);
       Monoalg_sim.original_discretization=mesh_resolution_hexa*10000;
       Monoalg_sim.desired_discretization=mesh_resolution_hexa*10000;

       meshfile=strcat('local_meshes/ieee_monoalg3d_2023/rdoste/',strcat('Fields_',name_final,'.alg'));

       Monoalg_ini_creation(name_mono,monodir,referenceMonoAlg,meshfile,MeshHex,Monoalg_sim,sigma,Electrodes_final_monoAlg,rootnodes)


       cd(strcat(resultspath))
end

