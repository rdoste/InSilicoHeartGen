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

function write_vtk_points(filename,point_coordinates,pointIDs,pointLabel)

%write points in vtk file for glyph visualization
fid = fopen(filename,'w','b');
fprintf(fid,'# vtk DataFile Version 3.0\x0A');
fprintf(fid,'vtk output\x0A');
fprintf(fid,'ASCII\x0A');
fprintf(fid,'DATASET UNSTRUCTURED_GRID\x0A');
str = sprintf('POINTS %d float\x0A',length(pointIDs));
fprintf(fid,str);
fprintf(fid,'%8.3f %8.3f %8.3f \n',(point_coordinates(pointIDs,:))');
fprintf(fid,'\n');

ex=exist('pointLabel');
if ex==1
fprintf(fid,'POINT_DATA %d\n',length(pointIDs));
fprintf(fid,'SCALARS Scalar1 double \n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,' %f \n',pointLabel);
end

fclose(fid);

end