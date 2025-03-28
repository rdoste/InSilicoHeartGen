function Electrodes_final=electrode_generation(Meshdata,referenceFolder)

%read ref data
Refdata0=load(strcat(referenceFolder,'\Reference.mat'));
Refdata=Refdata0.Reference;
%read electrode position
    fileID = fopen(strcat(referenceFolder,'\Electrodes.vtk'),'r');
    format='%f %f %f';
    Electrodes = cell2mat(textscan(fileID,format,10,'HeaderLines',5,'CollectOutput',1));
    fclose(fileID);



%obtain original coordinate frame
e1=normr(Refdata.a2p_vector);
e2=Refdata.r2l_vector;
e2=normr(e2-bsxfun(@times,dot(e1',e2')',e1));
e3=normr(cross(e1,e2));


%obtain original coordinate frame
n1=normr(Meshdata.a2p_vector);
n2=Meshdata.r2l_vector;
n2=normr(n2-bsxfun(@times,dot(n1',n2')',n1));
n3=normr(cross(n1,n2));

Re=normalize( [e1(:),e2(:),e3(:)] ,1,'norm');
Rn=normalize( [n1(:),n2(:),n3(:)],1,'norm');  
R=Re*Rn.';

%% Align heart
v2=(Refdata.v)*R;

 %heart center location in reference
septum_ref=find(abs(Refdata.r2l_geo)==0 );               
[k,dist] = dsearchn([Refdata.a2b_cobi(septum_ref),Refdata.r(septum_ref)],[0.5,5/6]);
center_point_Ref=v2(septum_ref(k),:);


%heart center location in mesh
septum_mesh=find(abs(Meshdata.r2l_geo)==0 );               
[k,dist] = dsearchn([Meshdata.a2b_cobi(septum_mesh),Meshdata.r(septum_mesh)],[0.5,5/6]);
center_point_Mesh=Meshdata.v(septum_mesh(k),:);

translate_vector=center_point_Mesh-center_point_Ref;

v2=v2+translate_vector;
% elem=Refdata.f;
% write_vtk_rbm('HeartAligned.vtk',length(v2),v2,length(elem),elem,ones(1,length(v2))); %write meshend


%% Align Electrodes

Electrodes_final=Electrodes*R+translate_vector;
write_vtk_points('Electrodes.vtk',Electrodes_final,1:10,1:10);




