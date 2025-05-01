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


function write_vtk_surf(varargin)
    %examples
    %write_vtk_surf('Torso_MLCC.vtk',pto_torso_dest,car_torso+1);
    %write_vtk_surf('Torso_MLCC.vtk',pto_torso_dest,car_torso+1,label);
    
    name=cell2mat(varargin(1));
    points_full=cell2mat(varargin(2));
    faces_full=cell2mat(varargin(3));
    N_points_full=length(points_full);
    N_faces_full=length(faces_full);
    
    
    fileID = fopen(name, 'w');
        fprintf(fileID,'# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET POLYDATA\n' );
        fprintf(fileID,'POINTS %d',N_points_full);
        fprintf(fileID,' float\n');
        fprintf(fileID,'%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',points_full');
        % fprintf(fileID,' \n');
        fprintf(fileID,'POLYGONS %d %d\n',N_faces_full, N_faces_full*4);
        faces_full2=[repmat(size(faces_full,2),size(faces_full,1),1) faces_full-1];
        fprintf(fileID,[repmat('%d\t', 1, size(faces_full2, 2)) '\n'], faces_full2');
        fprintf(fileID,' \n');
    if nargin==4
         label=cell2mat(varargin(4));
         if length(label)==length(faces_full)
             fprintf(fileID,'CELL_DATA %d\n',N_faces_full); 
             fprintf(fileID,'SCALARS labelcell double \n');
             fprintf(fileID,'LOOKUP_TABLE default\n');
             fprintf(fileID,[repmat('%f\t', 1, size(label, 2)) '\n'], label');
         elseif length(label)==length(points_full)
             fprintf(fileID,'POINT_DATA %d\n',N_points_full); 
             fprintf(fileID,'SCALARS labelpoint double 1\n');
             fprintf(fileID,'LOOKUP_TABLE default\n');
             fprintf(fileID,[repmat('%f\t', 1, size(points_full(:,1), 2)) '\n'],label');
         end
    end
     fclose(fileID);
 
end
