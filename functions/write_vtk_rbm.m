function write_vtk_rbm(varargin)
%function write_vtk_rbm(name,N_points,v,N_faces,f,label_1,label_2,label_3,fibers)
 
 
 
 name=cell2mat(varargin(1));
 N_points=cell2mat(varargin(2));
 v=cell2mat(varargin(3));
 N_faces=cell2mat(varargin(4));
 f=cell2mat(varargin(5));
 label_1=cell2mat(varargin(6));
 
 
 
 if nargin==6
     
     
     
      fileID = fopen(name,'w');
fprintf(fileID,'# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\n' );
fprintf(fileID,'POINTS %d',N_points);
fprintf(fileID,' float\n');
fprintf(fileID,'%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',v');
fprintf(fileID,'\n');
fprintf(fileID,'\n');
fprintf(fileID,'CELLS %d %d\n',N_faces,N_faces*5);
fprintf(fileID,'4 %d %d %d %d \n',(f(:,2:5))');
fprintf(fileID,'\n');
fprintf(fileID,'CELL_TYPES %d\n',N_faces);
fprintf(fileID,'%d \n',zeros(size(f(:,1)))+10);
fprintf(fileID,'\n');
fprintf(fileID,'POINT_DATA %d\n',N_points);
fprintf(fileID,'SCALARS Tpsi double \n');
fprintf(fileID,'LOOKUP_TABLE default\n');
fprintf(fileID,' %f \n',label_1);
fprintf(fileID,'\n');
fclose(fileID);  




elseif nargin ==7
 
 
 
 fibers=cell2mat(varargin(7));
 
     
     
          
     
fileID = fopen(name,'w');
fprintf(fileID,'# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\n' );
fprintf(fileID,'POINTS %d',N_points);
fprintf(fileID,' float\n');
fprintf(fileID,'%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',v');
fprintf(fileID,'\n');
fprintf(fileID,'\n');
fprintf(fileID,'CELLS %d %d\n',N_faces,N_faces*5);
fprintf(fileID,'4 %d %d %d %d \n',(f(:,2:5))');
fprintf(fileID,'\n');
fprintf(fileID,'CELL_TYPES %d\n',N_faces);
fprintf(fileID,'%d \n',zeros(size(f(:,1)))+10);
fprintf(fileID,'\n');
fprintf(fileID,'POINT_DATA %d\n',N_points);
fprintf(fileID,'SCALARS Scalar1 double \n');
fprintf(fileID,'LOOKUP_TABLE default\n');
fprintf(fileID,' %f \n',label_1);
fprintf(fileID,'\n');
fprintf(fileID,'FIELD FieldData 1 \n');
fprintf(fileID,'Fibers 3 %d double\n',N_points);
fprintf(fileID,'%f %f %f \n',fibers');
fclose(fileID);
     
     
  
     
 elseif nargin ==9
 
 
 label_2=cell2mat(varargin(7));
 label_3=cell2mat(varargin(8));
 fibers=cell2mat(varargin(9));
 
     
     
          
     
fileID = fopen(name,'w');
fprintf(fileID,'# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\n' );
fprintf(fileID,'POINTS %d',N_points);
fprintf(fileID,' float\n');
fprintf(fileID,'%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',v');
fprintf(fileID,'\n');
fprintf(fileID,'\n');
fprintf(fileID,'CELLS %d %d\n',N_faces,N_faces*5);
fprintf(fileID,'4 %d %d %d %d \n',(f(:,2:5))');
fprintf(fileID,'\n');
fprintf(fileID,'CELL_TYPES %d\n',N_faces);
fprintf(fileID,'%d \n',zeros(size(f(:,1)))+10);
fprintf(fileID,'\n');
fprintf(fileID,'POINT_DATA %d\n',N_points);
fprintf(fileID,'SCALARS Scalar1 double \n');
fprintf(fileID,'LOOKUP_TABLE default\n');
fprintf(fileID,' %f \n',label_1);
fprintf(fileID,'\n');
fprintf(fileID,'FIELD FieldData 3 \n');
fprintf(fileID,'Fibers 3 %d double\n',N_points);
fprintf(fileID,'%f %f %f \n',fibers');
fprintf(fileID,'Scalar2 1 %d',N_points);
fprintf(fileID,'double\n');
fprintf(fileID,'%f \n',label_2);
fprintf(fileID,'Scalar3 1 %d',N_points);
fprintf(fileID,'double\n');
fprintf(fileID,'%f \n',label_3);
fclose(fileID);
     
     
 elseif nargin ==8
 
 
 label_2=cell2mat(varargin(7));
 label_3=cell2mat(varargin(8));
 
     
     
fileID = fopen(name,'w');
fprintf(fileID,'# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\n' );
fprintf(fileID,'POINTS %d',N_points);
fprintf(fileID,' float\n');
fprintf(fileID,'%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',v');
fprintf(fileID,'\n');
fprintf(fileID,'\n');
fprintf(fileID,'CELLS %d %d\n',N_faces,N_faces*5);
fprintf(fileID,'4 %d %d %d %d \n',(f(:,2:5))');
fprintf(fileID,'\n');
fprintf(fileID,'CELL_TYPES %d\n',N_faces);
fprintf(fileID,'%d \n',zeros(size(f(:,1)))+10);
fprintf(fileID,'\n');
fprintf(fileID,'POINT_DATA %d\n',N_points);
fprintf(fileID,'FIELD FieldData 3 \n');
fprintf(fileID,'label1 1 %d double\n',N_points);
fprintf(fileID,'%f %f %f \n',label_1);
fprintf(fileID,'label2 1 %d',N_points);
fprintf(fileID,'double\n');
fprintf(fileID,'%f \n',label_2);
fprintf(fileID,'label3 1 %d',N_points);
fprintf(fileID,'double\n');
fprintf(fileID,'%f \n',label_3);
fclose(fileID);
 end
end
        