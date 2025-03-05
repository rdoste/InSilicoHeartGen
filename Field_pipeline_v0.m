clear
%% add folders with matlab libraries and functions
%set one drive path
current_path=pwd;
addpath (genpath(strcat(current_path,'\functions')));
addpath (genpath(strcat(current_path,'\dependencies')));

%% configure input data

origpath=strcat(current_path,'\inputs'); % path with the surface meshes
resultspath=strcat(current_path,'\outputs'); % path with the resulting files and fields
referenceFolder=strcat(current_path,'\PostDoc\UKBB\functions\activation_generic_cobi_CR'); %reference folder for activation and electrodes
referenceAlya=strcat(current_path,'\PostDoc\UKBB\functions\Alya_generic'); %reference folder for Alya files generation
referenceMonoAlg=strcat(current_path,'\PostDoc\MeshTool_paper\Rodero_Meshes\functions\Alya_generic'); %reference folder for ini files generation
        
Prefix_LV='LV_endo_ex_';   %change according to the format of the input data
Prefix_LV_epi='LV_epi_ex_';   %change according to the format of the input data
Prefix_RV='RV_ex_';   %change according to the format of the input data

meshformat='UKBB'  ; %type of mesh  UKBB--> inputs obtained from the UKBB image data
                   %               cut --> cut biventricular geomtry
                   %               biventricular_open --> biventricular  with open valves
                   %               biventricular_closed  --> biventricular  with closed valves (EM simulations) 
meshtype='all'    ; %type of mesh  all--> Tetrahedral (coarse and fine) and hexahedral
                   %               cut --> cut biventricular geomtry
                   %               biventricular_open --> biventricular  with open valves
mesh_resolution_fine=1;
mesh_resolution_coarse=1.5;
mesh_resolution_hexa=0.04;
                                 
%% get files name
cd(origpath)
cases=dir(strcat(Prefix_LV,'*'));
case_names=string({cases(:).name})';
case_numbers=str2double(extract(case_names,digitsPattern));
total_cases=length(case_numbers);
cd(resultspath)
global case_number ;

%% mesh generation
for index=1
    tStart = tic; 
        case_number=case_numbers(index);
         
        %% create case folder of input data
        if ~exist(num2str(case_number),'dir')
          mkdir(num2str(case_number));
          mkdir(num2str(case_number),'input');
          copyfile(strcat(origpath,'\',Prefix_LV,num2str(case_number),'.ply'),strcat(resultspath,'\',num2str(case_number),'\input'));
          copyfile(strcat(origpath,'\',Prefix_LV_epi,num2str(case_number),'.ply'),strcat(resultspath,'\',num2str(case_number),'\input'));
          copyfile(strcat(origpath,'\',Prefix_RV,num2str(case_number),'.ply'),strcat(resultspath,'\',num2str(case_number),'\input'));
        end
        
        
         %% read mesh

                cd(strcat(resultspath,'\',num2str(case_number),'\','input'))

                %1 create epicardium of right ventricle
                add_RV_width(Prefix_RV,3); 
                %2 fix LV mitral borders
                fix_LV_borders (Prefix_LV);
                [original_LV_mesh]=readVTK(strcat(Prefix_LV,num2str(case_number),'.ply'));%used for labelling

                %3 merge all surfaces
                labels0=merge_surfaces(Prefix_LV,Prefix_LV_epi,Prefix_RV);
                %4 add lid to LV
                add_lid_LV()
                %5 Remesh surfaces (important step, it fixes mesh issues)
                surf_edge_length=1.5;
                surf0=remesh_surfaces(surf_edge_length);

                %6 Move surfaces meshes to result folder and change directory to
                %that folder
                cd(strcat(resultspath,'\',num2str(case_number)))
                movefile(strcat(resultspath,'\',num2str(case_number),'\input\Closed_final_',num2str(case_number),'.ply'),strcat(resultspath,'\',num2str(case_number),'\'));
                movefile(strcat(resultspath,'\',num2str(case_number),'\input\labels0.vtk'),strcat(resultspath,'\',num2str(case_number),'\'));



                %% generate coarse mesh

                %find unit scale
                dist=sqrt(sum((surf0.points-repmat(surf0.points(1,:),length(surf0.points),1)).^2,2));
                scaledist=max(dist); %if ~100 --> mm

                if scaledist >50 &&  scaledist <500
                    disp('scale is in mm ')

                elseif scaledist >5 &&  scaledist <=50
                    disp('scale is in cm -->Units for Alya and Personalization ')
                    disp('conversion to mm for meshing')
                     surf0.points=surf0.points.*10;
                     labels0.points=labels0.points.*10;

                else
                    error('check your mesh dimensions (not mm nor cm)')
                end
         % if ~exist('coarse.vtu','file')    
                %coarse tetrahedral  mesh at 1.5 mm resolution
                MeshCoarse=tetrahedral_meshing_cgal('labels0.vtk',surf0,1/4,mesh_resolution_coarse);

            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            disp('conversion to cm after meshing -->Units for Alya and Personalization ')
            MeshCoarse.points=MeshCoarse.points./10;
            surf0.points=surf0.points./10;
            disp('units converted to cm')
            vtkWrite(MeshCoarse, 'Coarse.vtu');



        % else
        %    MeshCoarse=vtkRead('Coarse.vtu');
        % 
        % end

        %% generate fine mesh at 1 mm resolution
        disp('generating fine mesh with mmg')
         % if ~exist('fine.vtu','file')
            Meshmm=MeshCoarse;
            Meshmm.points=MeshCoarse.points.*10;
            mmgSizingParam = [0.1 0.9 1.1];
            meanEdgLen = mesh_resolution_fine;
            mmgSizingParam2 =[mmgSizingParam(1) mmgSizingParam(2:3) *meanEdgLen];
            [MeshFine,mmgStatus,mmgOutput] = mmg_RD(Meshmm, sprintf(' -hausd %1.5e -hmin %1.5e -hmax %1.5e',  mmgSizingParam2(:)'));
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            disp('conversion to cm after meshing -->Units for Alya and Personalization ')
            MeshFine.points=MeshFine.points./10;
            disp('units converted to cm')
            vtkWrite(MeshFine, 'fine.vtu');

         % else
         %    MeshFine=vtkRead('fine.vtu');
         % end




         %% generate Hexa mesh at 1 mm resolution
          
        [ALG,tet_ID,bar]= hexa_mesher('Coarse.vtu', mesh_resolution_hexa,1);
        outputmesh=movefile('hex_Coarse.vtk',strcat('Hexa_',num2str(case_number),'.vtk'));
        MeshHex=vtkRead(strcat('Hexa_',num2str(case_number),'.vtk'));



        
        % generate labels
         % if ~exist('labels_final.vtk','file')
            disp('generating labels')
            cd(strcat(resultspath,'\',num2str(case_number),'\input\'))        
            sur_coarse = vtkDataSetSurfaceFilter(MeshCoarse);
            original_LV_mesh.points=original_LV_mesh.points./10;

            labelfinal3=Labels_UKBB(sur_coarse,original_LV_mesh);

            movefile(strcat(resultspath,'\',num2str(case_number),'\input\labels_final.vtk'),strcat(resultspath,'\',num2str(case_number),'\'));
         % end
       

        
        
        % generate fields
        cd(strcat(resultspath,'\',num2str(case_number)))

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
        Field_generator_UKBB_function100(Fiber_info,meshformat,pericardium_level, epiendo, epiendoRV);
        % toc

        cd(strcat(resultspath,'\',num2str(case_number)))

        %% creation of files from Alya
        nameAlyaFolder=strcat('Alya_Control_',num2str(case_number));


        Alyafilespath=strcat(resultspath,'\',num2str(case_number),'\',nameAlyaFolder);
        Data4Alya=strcat(resultspath,'\',num2str(case_number),'\ensi_Fine_',num2str(case_number),'\Case_Fine.mat') ; %choose between the mesh files  (coarse for EM or fine for EP)
        Electrodes_final=electrode_generation(Data4Alya,referenceFolder);

        name_alya=strcat('heart_',num2str(case_number));
        AlyaFiles_creation_all(Alyafilespath,Data4Alya,name_alya,referenceAlya,Electrodes_final,case_number);

        %%generate activation
        BCL=0.8;
        generateActivationAlya_Closed_v2(name_alya,Data4Alya,referenceFolder,Alyafilespath,'generic',BCL);



%         %% Hexa Files generation
%         %read field data
        casepath=strcat(resultspath,'/',num2str(case_number),'\',strcat('ensi',num2str(case_number)));
        monodir=strcat(resultspath,'\',num2str(case_number),'\MonoAlg3D');

        cd(casepath)
        Data=load('Case_coarse.mat');
        %read mesh data

       cd(strcat(resultspath,'\',num2str(case_number)))

        mkdir(monodir)
        cd(monodir);
        %generate ALG file
        HexaFieldsGeneration_function_cells_v2(monodir,ALG,MeshCoarse,MeshHex,Data,tet_ID,bar,epiendo,epiendoRV);


        %new electrodes
        Electrodes_orig=electrode_generation_v2(Data,referenceFolder);
        Electrodes_final=[Electrodes_orig(:,1)-min(Data.v(:,1)),Electrodes_orig(:,2)-min(Data.v(:,2)),Electrodes_orig(:,3)-min(Data.v(:,3))];
        write_vtk_points('Electrodes_final.vtk',Electrodes_final,1:10,1:10);

       %root nodes
       roots_number=4;
       BCL=0.8;
       name_mono=strcat('UKBB_',num2str(case_number));
       rootnodes=rootnodes_from_IDs(name_mono,Data,referenceFolder,monodir,roots_number,BCL);
       rootnodes.time=rootnodes.time-0.02;  %transformation from Alya generic (simulations starts at 0.02)
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

       meshfile=strcat('local_meshes/ieee_monoalg3d_2023/rdoste/',strcat('Fields_',num2str(case_number),'.alg'));

       Monoalg_ini_creation(name_mono,monodir,referenceMonoAlg,meshfile,MeshHex,Monoalg_sim,sigma,Electrodes_final,rootnodes)


       cd(strcat(resultspath))
end
    tfinal=toc(tabs);
    %sprintf('Total time is %s ',tfinal)
