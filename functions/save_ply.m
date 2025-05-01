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

function save_ply(pto,car,name)

    nnode=length(pto);
    nelem=length(car);
    
    fileID = fopen(strcat(name,'.ply'),'w');
    fprintf(fileID,'ply\nformat ascii 1.0\ncomment VTK generated PLY File\nobj_info vtkPolyData points and polygons: vtk4.0\n' );
    fprintf(fileID,'element vertex %d',nnode);
    fprintf(fileID,'\nproperty float x\nproperty float y\nproperty float z\n' );
    fprintf(fileID,'element face %d',nelem);
    fprintf(fileID,'\nproperty list uchar int vertex_indices\nend_header\n' );
    fprintf(fileID,'%.4f %.4f %.4f \n',pto');
    fprintf(fileID,'3 %d %d %d \n',car');
    fclose(fileID);

end