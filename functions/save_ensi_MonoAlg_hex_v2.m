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

function save_ensi_MonoAlg_hex_v2(v,f,Ventricle,d3,Tphi3,Tphi_cobi,Epiendo,Epiendo3,a2b_uvc,a2b_cobi,r,lvrv_cobi,r2l_geo,a2b,r2l,a2p,F,F_S,F_N,Material,aha,a2b_projected,tm_projected,tv_projected,rt_projected)
 %function save_ensi_MonoAlg(v,f,Ventricle,d3)

%create ensight file
%create geo
save_ensi_geo_hex('Fields_Hex',v,f);

%create fields
save_ensi_field_hex('Fields_Hex','Ventricle',Ventricle);
save_ensi_field_hex('Fields_Hex','Transmurality',d3);
save_ensi_field_hex('Fields_Hex','TransmuralityRV',Tphi3);
save_ensi_field_hex('Fields_Hex','Transmurality_cobi',Tphi_cobi);
save_ensi_field_hex('Fields_Hex','celltype',Epiendo);
save_ensi_field_hex('Fields_Hex','celltypeRVsept',Epiendo3);
save_ensi_field_hex('Fields_Hex','apex_2_base_uvc',a2b_uvc);
save_ensi_field_hex('Fields_Hex','apex_2_base_cobi',a2b_cobi);
save_ensi_field_hex('Fields_Hex','Circunferential_cobi',r);
save_ensi_field_hex('Fields_Hex','LVvsRV_cobi',lvrv_cobi);
save_ensi_field_hex('Fields_Hex','rv_2_lv_geodesic',r2l_geo);
save_ensi_field_hex('Fields_Hex','apex_2_base_projection',a2b);
save_ensi_field_hex('Fields_Hex','rv_2_lv_projection',r2l);
save_ensi_field_hex('Fields_Hex','anterior_2_posterior_projection',a2p);
save_ensi_field_hex('Fields_Hex','FiberL',F);
save_ensi_field_hex('Fields_Hex','FiberS',F_S);
save_ensi_field_hex('Fields_Hex','FiberN',F_N);
save_ensi_field_hex('Fields_Hex','Material',Material);
save_ensi_field_hex('Fields_Hex','aha',aha);

% %create case
fields_cell{1}={'Ventricle','Ventricle',1}; %name of variable in matlab, name for variable in ensigth, dimension
fields_cell{2}={'Transmurality','Transmurality',1}; %name of file, name for variable in ensigth, dimension
fields_cell{3}={'TransmuralityRV','TransmuralityRV',1}; %name of file, name for variable in ensigth, dimension
fields_cell{4}={'Transmurality_cobi','Cobi_Transmurality',1}; %name of file, name for variable in ensigth, dimension
fields_cell{5}={'celltype','celltype',1}; %name of file, name for variable in ensigth, dimension
fields_cell{6}={'celltypeRVsept','celltypeRVsept',1}; %name of file, name for variable in ensigth, dimension
fields_cell{7}={'apex_2_base_uvc','apex_2_base_uvc',1}; %name of file, name for variable in ensigth, dimension
fields_cell{8}={'apex_2_base_cobi','Cobi_apex_2_base',1}; %name of file, name for variable in ensigth, dimension
fields_cell{9}={'Circunferential_cobi','Cobi_Circunferential',1}; %name of file, name for variable in ensigth, dimension
fields_cell{10}={'LVvsRV_cobi','Cobi_LVvsRV',1}; %name of file, name for variable in ensigth, dimension
fields_cell{11}={'rv_2_lv_geodesic','rv_2_lv_geodesic',1}; %name of file, name for variable in ensigth, dimension
fields_cell{12}={'apex_2_base_projection','apex_2_base_projection',1}; %name of file, name for variable in ensigth, dimension
fields_cell{13}={'rv_2_lv_projection','rv_2_lv_projection',1}; %name of file, name for variable in ensigth, dimension
fields_cell{14}={'anterior_2_posterior_projection','anterior_2_posterior_projection',1}; %name of file, name for variable in ensigth, dimension
fields_cell{15}={'FiberL','FiberL',3}; %name of file, name for variable in ensigth, dimension
fields_cell{16}={'FiberS','FiberS',3}; %name of file, name for variable in ensigth, dimension
fields_cell{17}={'FiberN','FiberN',3}; %name of file, name for variable in ensigth, dimension
fields_cell{18}={'Material','Material',1}; %name of file, name for variable in ensigth, dimension
fields_cell{19}={'aha','aha',1}; %name of file, name for variable in ensigth, dimension



save_ensi_case_hex('Fields_Hex',fields_cell);


%functions
    function save_ensi_geo_hex(name,node_c,elem_c)
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
        fwrite(fid,                 N_elem_c, '*uint32','l');
        fwrite(fid,(1:length(elem_c)),'*uint32','l' );
        fwrite(fid,elem_c_write+1,'*uint32','l' );
        fprintf(fid,'%-80s',             'tetra4' );
        fwrite(fid,                 0, '*uint32','l');
        fclose(fid);
    
    
    end

    function save_ensi_field_hex(name,namefield,Field)
      % save_ensi_field('Test','phi',phi)

        filename=strcat(name,'.ensi.',namefield);     
        fid=fopen(filename,'w');
        fprintf(fid,'%-80s','description line 1');
        fprintf(fid,'%-80s','part');
        fwrite(fid,                1, '*uint32','l');
        fprintf(fid,'%-80s',            'hexa8');
        fwrite(fid,Field(:),'*float','l' );
        fclose(fid);
    
    end

    function save_ensi_case_hex(name,fieldscell)
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
            if  fieldscell{1,ind}{1,3}==1
            fprintf(fid,'scalar per element: ');
            elseif fieldscell{1,ind}{1,3}==3
             fprintf(fid,'vector per element: ');
            end
            fprintf(fid,'%d',fieldscell{1,ind}{1,3});
            fprintf(fid,' %s',fieldscell{1,ind}{1,2});
            fprintf(fid,' %s \n',strcat(Fieldname,fieldscell{1,ind}{1,1}));
            
        end
        fclose(fid);
    end

end