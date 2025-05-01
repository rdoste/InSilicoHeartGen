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
%     along with this program.  If not, see <https://www.gnu.org/licenses/>

function   [N_points_full,pto,N_faces_full,car,Fid]=read_labels2(name)

    fileID = fopen(name,'r');
   %Read number of points
    format='%s %f %s';
    data = textscan(fileID,format,1,'HeaderLines',4,'Delimiter',' ');
    N_points_full=data{1,2};
   %save points
    format=' %f %f %f';
    data= textscan(fileID,format,N_points_full,'CollectOutput',1);
    pto=data{1,1};
   %Read number of faces
    format='%s %f %f';
    data = textscan(fileID,format,1,'Delimiter',' ');
    
    N_faces_full=data{1,2}(1);
    if isnan(N_faces_full)==1
       data = textscan(fileID,format,1,'Delimiter',' ');
    end
        
    connectivity=data{1,3}/N_faces_full;
   %save faces
    format=repmat('%d ',1,connectivity);
    data= textscan(fileID,format,N_faces_full,'CollectOutput',1);
    car=data{1}(:,2:end);
    
    
    format='%f';
    data= textscan(fileID,format,N_faces_full,'HeaderLines',4);
    k=1;
    while isempty(data{1,1})==1 && k<50
       data= textscan(fileID,format,N_faces_full,'HeaderLines',k);
       k=k+1;
    end
        
    Fid=data{1,1};
    
    fclose(fileID);

end