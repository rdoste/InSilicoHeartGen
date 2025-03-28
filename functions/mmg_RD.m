function [mesh,status,cmdout] = mmg_RD(mesh,  paramString)

[tmpdir,name] = fileparts(tempname);
meshfile = sprintf('%s/%s.mesh', tmpdir, name);
solfile = sprintf('%s/%s.sol', tmpdir, name);

mmgWriteMesh(mesh, meshfile);

mpath = "C:\Users\rubste\Documents\Cobiveco\functions";
  if isfolder('C:\DATA\')% personal laptop
    mpath="C:\DATA\PostDoc\Cobiveco\functions";
  end

[status,cmdout] = system(sprintf("%s/../dependencies/mmg/build/bin/mmg3d_O3 %s %s", mpath, string(meshfile), string(paramString)));
cmdout
mesh = struct();
 if status==0
        mesh = mmgReadMesh(strcat(meshfile(1:end-5),'.o.mesh'));
 end

end
