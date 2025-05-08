%     Automated pipeline for large-scale cardiac in silico trials 
%     Copyright (C) 2025 Ruben Doste. Contact: ruben.doste@gmail.com
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

function [interp, interp2, interp3]= calculate_ray_intersections(points_full,faces_full,vertex1,vertex2,vertex3)

%computes the number of ray intersectiosn from faces normal


%% Parameters;

        %%

         TR=triangulation(faces_full,points_full);
         facenormal0 = faceNormal(TR);
         centroid=meshcentroid(points_full,faces_full);
         length_target=length(vertex1);
          

        %% point ray tracing

            %USING FACES
                % pointnormal=vertexNormal(TR);
                % pointnormal=gpuArray(pointnormal);
                

         disp('performing ray tracing algorithm...');

                interp=zeros(size(centroid,1),1);
                interp2=interp;
                interp3=interp;

                tic
             checkGPU=gpuDeviceCount;
             if checkGPU~=0
                 facenormal = gpuArray(facenormal0);
                 % vertex1=vertex1';
                 % vertex2=vertex2';
                 % vertex3=vertex3';

                 for indxp=1:length(centroid)
                    [intersection,~]= arrayfun(@rayTriGPU,vertex1(:,1)',vertex1(:,2)',vertex1(:,3)',...
                        vertex2(:,1)',vertex2(:,2)',vertex2(:,3)',...
                        vertex3(:,1)',vertex3(:,2)',vertex3(:,3)',...
                        centroid(indxp,1),centroid(indxp,2),centroid(indxp,3),...
                        facenormal(indxp,1),facenormal(indxp,2),facenormal(indxp,3));


                    interp3(indxp) = sum(gather(intersection)>0);
                    interp(indxp) = min(intersection);
                    intersection(intersection<-1e-10)=999;
                    interp2(indxp) = abs(min(intersection));


                 end

             else

                parfor ind_origin=1:length(centroid)
                    [TTT1,intersection1]= TriangleRayIntersection ( repmat(centroid(ind_origin,:),length_target,1), repmat(facenormal0(ind_origin,:),length_target,1), vertex1, vertex2, vertex3,'planeType','two sided');
                    [TTT2]= TriangleRayIntersection ( repmat(centroid(ind_origin,:),length_target,1), repmat(facenormal0(ind_origin,:)*-1,length_target,1), vertex1, vertex2, vertex3,'planeType','two sided');

                    interp3(ind_origin) = sum(intersection1(or(TTT2,TTT1))>0);
                    interp(ind_origin) = min(intersection1(or(TTT2,TTT1)));
                    intersection1(intersection1<-1e-10)=999;
                    interp2(ind_origin) = abs(min(intersection1(or(TTT2,TTT1))));

                end

                
             end
                interp=interp';
                interp2=interp2';
                interp3=interp3';
                 toc

                 

 

end
