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


function     Gradient=gradientinterpol(Gradient_ori,node_iris,pointsTetra,bar)
             %Field interpolation using baricentric coordinates and GPU
             
             checkGPU=gpuDeviceCount;
             if checkGPU~=0
                        DT=Gradient_ori(pointsTetra,:);
                        DT1=reshape(DT',3,node_iris,4);
                        DT2=permute(DT1,[3 1 2]);
                        M1= gpuArray(DT2);
                        bar2=reshape(bar',1,4,node_iris);
                        bar22=repmat(bar2,3,1,1);
                        bar222=permute( bar22,[2 1 3]);
                        M2= gpuArray(bar222);
                        Q2= pagefun(@times, M1, M2);
                        Q2= gather(Q2);
                        Q3=sum(Q2,1);
                        Gradient=reshape(Q3,3,node_iris)';
                        
                        gpuDevice(1);
             else
                        DT=Gradient_ori(pointsTetra,:);
                        DT1=reshape(DT',3,node_iris,4);
                        DT2=permute(DT1,[3 1 2]);                     
                        bar2=reshape(bar',1,4,node_iris);
                        bar22=repmat(bar2,3,1,1);
                        bar222=permute( bar22,[2 1 3]);                              
                             
                        QX2=zeros(size(DT2));        
                        for indx=1:length(bar)
                               QX2(:,:,indx)= times(DT2(:,:,indx), bar222(:,:,indx));
                        end
                        QX3=sum(QX2,1);        
                        GR=QX3(1,:,:);
                        Gradient=reshape(GR,3,node_iris)';
             end


end