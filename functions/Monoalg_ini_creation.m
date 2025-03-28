function Monoalg_ini_creation(name,monodir,referenceMonoAlg,meshfile,MeshHex,Monoalg_sim,sigma,Electrodes_final,rootnodes)%Alya Files creatio

%automatic MonoAlg3D files generation

 
    
    %% generation of ini file

%%.ini-------------------------------------------------------------------------------------------
copyfile(strcat(referenceMonoAlg,'\','MonoAlg3D_generic.ini'),strcat(monodir,'\',name,'.ini'));
fid  = fopen(strcat(name,'.ini'),'r');
f=fread(fid,'*char')';
fclose(fid);
f = strrep(f,'**num_threads**',num2str(Monoalg_sim.num_threads));
f = strrep(f,'**dt_pde**',num2str(Monoalg_sim.dt_pde));
f = strrep(f,'**duration**',num2str(Monoalg_sim.duration));
f = strrep(f,'**output_dir**',Monoalg_sim.output_dir);
f = strrep(f,'**sigma_l**',num2str(sigma(1)));
f = strrep(f,'**sigma_t**',num2str(sigma(2)));
f = strrep(f,'**sigma_n**',num2str(sigma(3)));
f = strrep(f,'**fast_endo_layer_scale**',num2str(Monoalg_sim.fast_endo_layer_scale));
f = strrep(f,'**dt_ode**',num2str(Monoalg_sim.dt_ode));
f = strrep(f,'**meshfile**',meshfile);
f = strrep(f,'**num_volumes**',num2str(Monoalg_sim.num_volumes));
f = strrep(f,'**original_discretization**',num2str(Monoalg_sim.original_discretization));
f = strrep(f,'**desired_discretization**',num2str(Monoalg_sim.desired_discretization));


fid  = fopen(strcat(name,'.ini'),'w');
fprintf(fid,'%s',f);
fclose(fid);

%% add stimulus

for ind=1:length(rootnodes.IDs)
    fid  = fopen(strcat(name,'.ini'),'a');
    fprintf(fid,'%s\n',strcat( '[stim_root_node_',num2str(ind),']'));
    fprintf(fid,'start = %.2f \n',(rootnodes.time(ind)*1000));
    fprintf(fid,'duration = %.1f \n',4.0);
    fprintf(fid,'period = %.1f \n',Monoalg_sim.BCL);
    fprintf(fid,'current = %s \n',num2str(-53));
    fprintf(fid,'center_x = %.1f \n',rootnodes.coords(ind,1)*10000);
    fprintf(fid,'center_y = %.1f \n',rootnodes.coords(ind,2)*10000);
    fprintf(fid,'center_z = %.1f \n',rootnodes.coords(ind,3)*10000);
    fprintf(fid,'radius = %s \n',num2str(2000.0));
    fprintf(fid,'main_function = stim_sphere \n');
    fprintf(fid,'\n');
    
    fclose(fid);

end


%% add electrodes
     fid  = fopen(strcat(name,'.ini'),'a');
     fprintf(fid,'[calc_ecg] \n');
     fprintf(fid,'main_function=pseudo_bidomain \n');
     fprintf(fid,'init_function=init_pseudo_bidomain \n');   
     fprintf(fid,'end_function=end_pseudo_bidomain \n');
     fprintf(fid,'calc_rate=%d \n',length(Electrodes_final));
     for ind=1:length((Electrodes_final))
        fprintf(fid,'%s%.2f,%.2f,%.2f  \n',strcat('lead',num2str(ind),'='),Electrodes_final(ind,1)*10000,Electrodes_final(ind,2)*10000,Electrodes_final(ind,3)*10000);
     end
     fprintf(fid,'sigma_b = %s \n',num2str(20));
     fprintf(fid,'use_gpu=true \n');
     fclose(fid);





end
