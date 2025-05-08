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

function Monoalg_ini_creation(name,monodir,referenceMonoAlg,meshfile,MeshHex,Monoalg_sim,sigma,Electrodes_final,rootnodes)%Monoalg Files creatio

%automatic MonoAlg3D files generation

 
    
    %% generation of ini file

%%.ini-------------------------------------------------------------------------------------------
copyfile(fullfile(referenceMonoAlg,'MonoAlg3D_generic.ini'),fullfile(monodir,[name,'.ini']));
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
