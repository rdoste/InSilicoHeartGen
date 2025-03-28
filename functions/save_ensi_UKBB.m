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
function save_ensi_UKBB(v,f,Ventricle,d3,Tphi3,Tphi_cobi,Epiendo,Epiendo3,a2b_uvc,a2b_cobi,r,lvrv_cobi,r2l_geo,a2b,r2l,a2p,F,F_S,F_N,apex_2_base,aha,Material)

%function that saves the required fields in ensi format

%create ensight file
%create geo
save_ensi_geo('Fields',v,f);

%create fields
save_ensi_field('Fields','Ventricle',Ventricle);
save_ensi_field('Fields','Transmurality',d3);
save_ensi_field('Fields','TransmuralityRV',Tphi3);
save_ensi_field('Fields','Transmurality_cobi',Tphi_cobi);
save_ensi_field('Fields','celltype',Epiendo);
save_ensi_field('Fields','celltypeRVsept',Epiendo3);
save_ensi_field('Fields','apex_2_base_uvc',a2b_uvc);
save_ensi_field('Fields','apex_2_base_cobi',a2b_cobi);
save_ensi_field('Fields','Circunferential_cobi',r);
save_ensi_field('Fields','LVvsRV_cobi',lvrv_cobi);
save_ensi_field('Fields','rv_2_lv_geodesic',r2l_geo);
save_ensi_field('Fields','apex_2_base_projection',a2b);
save_ensi_field('Fields','rv_2_lv_projection',r2l);
save_ensi_field('Fields','anterior_2_posterior_projection',a2p);
save_ensi_field('Fields','FiberL',F);
save_ensi_field('Fields','FiberS',F_S);
save_ensi_field('Fields','FiberN',F_N);
save_ensi_field('Fields','Repolarisation_gradients',apex_2_base);
save_ensi_field('Fields','aha',aha);
save_ensi_field_tet('Fields','Material',Material);


%create case
fields_cell{1}={'Ventricle','Ventricle',size(Ventricle,2),size(Ventricle,1)}; %name of variable in matlab, name for variable in ensigth, dimension
fields_cell{2}={'Transmurality','Transmurality',size(d3,2),size(d3,1)}; %name of file, name for variable in ensigth, dimension
fields_cell{3}={'TransmuralityRV','TransmuralityRV',size(Tphi3,2),size(Tphi3,1)}; %name of file, name for variable in ensigth, dimension
fields_cell{4}={'Transmurality_cobi','Cobi_Transmurality',size(Tphi_cobi,2),size(Tphi_cobi,1)}; %name of file, name for variable in ensigth, dimension
fields_cell{5}={'Celltype','Celltype',size(Epiendo,2),size(Epiendo,1)}; %name of file, name for variable in ensigth, dimension
fields_cell{6}={'CelltypeRVsept','CelltypeRVsept',size(Epiendo3,2),size(Epiendo3,1)}; %name of file, name for variable in ensigth, dimension
fields_cell{7}={'Apex_2_base_uvc','Apex_2_base_uvc',size(a2b_uvc,2),size(a2b_uvc,1)}; %name of file, name for variable in ensigth, dimension
fields_cell{8}={'Apex_2_base_cobi','Apex_2_base_cobi',size(a2b_cobi,2),size(a2b_cobi,1)}; %name of file, name for variable in ensigth, dimension
fields_cell{9}={'Circunferential_cobi','Cobi_Circunferential',size(r,2),size(r,1)}; %name of file, name for variable in ensigth, dimension
fields_cell{10}={'LVvsRV_cobi','Cobi_LVvsRV',size(lvrv_cobi,2),size(lvrv_cobi,1)}; %name of file, name for variable in ensigth, dimension
fields_cell{11}={'Rv_2_lv_geodesic','Rv_2_lv_geodesic',size(r2l_geo,2),size(r2l_geo,1)}; %name of file, name for variable in ensigth, dimension
fields_cell{12}={'Apex_2_base_projection','Apex_2_base_projection',size(a2b,2),size(a2b,1)}; %name of file, name for variable in ensigth, dimension
fields_cell{13}={'Rv_2_lv_projection','Rv_2_lv_projection',size(r2l,2),size(r2l,1)}; %name of file, name for variable in ensigth, dimension
fields_cell{14}={'Anterior_2_posterior_projection','Anterior_2_posterior_projection',size(a2p,2),size(a2p,1)}; %name of file, name for variable in ensigth, dimension
fields_cell{15}={'FiberL','FiberL',size(F,2),size(F,1)}; %name of file, name for variable in ensigth, dimension
fields_cell{16}={'FiberS','FiberS',size(F_S,2),size(F_S,1)}; %name of file, name for variable in ensigth, dimension
fields_cell{17}={'FiberN','FiberN',size(F_N,2),size(F_N,1)}; %name of file, name for variable in ensigth, dimension
fields_cell{18}={'Repolarisation_gradients','Repolarisation_gradients',size(apex_2_base,2),size(apex_2_base,1)}; %name of file, name for variable in ensigth, dimension
fields_cell{19}={'Aha','Aha',size(aha,2),size(aha,1)}; %name of file, name for variable in ensigth, dimension
fields_cell{20}={'Material','Material',size(Material,2),size(Material,1)}; %name of file, name for variable in ensigth, dimension



save_ensi_case('Fields',fields_cell,length(v),length(f));


%% define functions
    % 
    function save_ensi_geo(name,node_c,elem_c)
      % save_ensi_geo('Test',v,f)
        
        elem_c_write=elem_c(:,2:end)';
        elem_c_write=elem_c_write(:);
        N_elem_c=length(elem_c);        
        N_nodes=length(node_c);

        namegeo=strcat(name,'.ensi.geo');
        fid=fopen(namegeo,'w');
        fprintf(fid,'%-80s','C Binary');
        fprintf(fid,'%-80s','description line 1');
        fprintf(fid,'%-80s','description line 2');
        fprintf(fid,'%-80s','node id given');
        fprintf(fid,'%-80s','element id given');
        fprintf(fid,'%-80s','part');
        fwrite(fid,                 1, '*uint32','l');
        fprintf(fid,'%-80s',            'description line 1' );
        fprintf(fid,'%-80s',            'coordinates');
        fwrite(fid,                 N_nodes, '*uint32','l');
        fwrite(fid,(1:length(node_c)),'uint32' );
        fwrite(fid,node_c(:,1),'float','l' );
        fwrite(fid,node_c(:,2),'float','l' );
        fwrite(fid,node_c(:,3),'float','l' );
        fprintf(fid,'%-80s',            'hexa8');
        fwrite(fid,                 0, '*uint32','l');
        fprintf(fid,'%-80s',             'tetra4' );
        fwrite(fid,                 N_elem_c, '*uint32','l');
        fwrite(fid,(1:length(elem_c)),'*uint32','l' );
        fwrite(fid,elem_c_write+1,'*uint32','l' );
        fclose(fid);
    
    end
    
    function save_ensi_field(name,namefield,Field)
    % save_ensi_field('Test','phi',phi)
        
        filename=strcat(name,'.ensi.',namefield);    
        fid=fopen(filename,'w');
        fprintf(fid,'%-80s','description line 1');
        fprintf(fid,'%-80s','part');
        fwrite(fid,                1, '*uint32','l');
        fprintf(fid,'%-80s',            'coordinates');
        fwrite(fid,Field(:),'*float','l' );
        fclose(fid);
    
    end

    function save_ensi_field_tet(name,namefield,Field)
        % save_ensi_field('Test','phi',phi)    
    
        filename=strcat(name,'.ensi.',namefield);       
        fid=fopen(filename,'w');
        fprintf(fid,'%-80s','description line 1');
        fprintf(fid,'%-80s','part');
        fwrite(fid,                1, '*uint32','l');
        fprintf(fid,'%-80s',            'tetra4');
        fwrite(fid,Field(:),'*float','l' );
        fclose(fid);
    
    end

    function save_ensi_case(name,fieldscell,length_v,length_f)    
         % save_ensi_case('Test','phi',phi)
    
        filename=strcat(name,'.ensi.case');
        geoname=strcat(name,'.ensi.geo');
        Fieldname=strcat(name,'.ensi.');
        
        fid=fopen(filename,'w');
        fprintf(fid,'FORMAT \n');
        fprintf(fid,'type:	ensight gold \n');
        fprintf(fid,'GEOMETRY \n');
        fprintf(fid,'model: 1 ');
        fprintf(fid,' %s \n',geoname);
        fprintf(fid,'VARIABLE \n');
        for ind=1:length(fieldscell)
            if  fieldscell{1,ind}{1,3}==1 && fieldscell{1,ind}{1,4}==length_v
            fprintf(fid,'scalar per node: ');
            elseif fieldscell{1,ind}{1,3}==3 && fieldscell{1,ind}{1,4}==length_v
             fprintf(fid,'vector per node: ');
            elseif fieldscell{1,ind}{1,3}==1 && fieldscell{1,ind}{1,4}==length_f
             fprintf(fid,'scalar per element: ');
            elseif fieldscell{1,ind}{1,3}==3 && fieldscell{1,ind}{1,4}==length_f
             fprintf(fid,'vector per element: ');    
            end
            fprintf(fid,'%d',fieldscell{1,ind}{1,3});
            fprintf(fid,' %s',fieldscell{1,ind}{1,2});
            fprintf(fid,' %s \n',strcat(Fieldname,fieldscell{1,ind}{1,1}));            
        end
        fclose(fid);
    end


end