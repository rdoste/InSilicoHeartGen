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

function fix_LV_borders (Prefix_LV, case_number)
    %read LV
        [pto,faces]=read_ply(strcat(Prefix_LV,num2str(case_number),'.ply'));   %surface mesh with no topological defects (can be the same as labels)
        faces_mod=faces;
        ind=1;
        while ind<3  
            boundary=compute_boundary(faces_mod);
            %remove boundary points and faces
            [indx_faces]=~ismember(faces_mod,boundary);
            faces_mod=faces_mod.*indx_faces;
            faces_mod= faces_mod(all(faces_mod,2),:);
        
            ind=ind+1;
        end 
    
    %save 
     save_ply(pto,faces_mod-1,strcat(Prefix_LV,num2str(case_number),'_mod'));

end


