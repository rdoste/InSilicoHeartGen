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